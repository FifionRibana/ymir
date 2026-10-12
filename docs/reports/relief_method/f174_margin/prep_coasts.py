"""ADR Finding 174-A2 -- the reference coast windows from Copernicus DEM GLO-30, as F161's Corsica grids.

Usage: python -I prep_coasts.py <tile dir> <out dir>

For each window (declared in `f174_declared.md`): the 1-arcsecond tile's pixel centres (PixelIsPoint: pixel (i, j) at
lon0 + i/w, lat0 + 1 - j/h) are projected to a local metric frame centred on the window (metres per degree from the WGS84
ellipsoid at each pixel's latitude) and averaged into a 1024 x 1024 grid of 48.83 m cells (50 km), rows SOUTH-FIRST: a cell
with no land pixel (DEM <= 0) is sea at -50 m, land keeps the mean of its land pixels -- F161's `prep_corse.py` exactly. The
fine grid is block-averaged twice to 256 x 256 cells of 195.3125 m, as `corse_1024` (sea and land averaged together). Every
island is kept (`prep_corse.py` dropped all but Corsica). The grid is written as <name>_256.bin (u32 w, u32 h, f32 LE).

Corsica west is a 256 x 256 crop of `data/corsica/corse_1024.bin` (the same 195.3125 m cells, frame centred at 42.18 N
9.045 E), centred on 42.20 N 8.65 E. A window is a mosaic of every 1-degree tile its 50 km box touches; a tile
Copernicus does not publish (HTTP 404: no land in it) is open sea.
"""
import sys
from pathlib import Path

import numpy as np
import tifffile

# name, type, centre latitude, centre longitude
WINDOWS = [
    ("galicia_rias_baixas", "rocky", 42.40, -8.70),
    ("dalmatia_split", "rocky", 43.55, 16.35),
    ("bretagne_nord", "rocky", 48.72, -3.65),
    ("landes_mimizan", "low", 44.20, -1.35),
    ("vendee_sables", "low", 46.55, -1.65),
    ("jutland_west", "low", 56.05, 8.30),
    ("golfe_du_lion_agde", "low", 43.25, 3.45),
]
N, NF = 256, 1024
CELL = 50_000.0 / NF
HALF = NF * CELL / 2


def m_per_deg(lat):
    p = np.radians(lat)
    return 111132.954 - 559.822 * np.cos(2 * p) + 1.175 * np.cos(4 * p), 111412.84 * np.cos(p) - 93.5 * np.cos(3 * p)


def tile_name(lat0, lon0):
    ns = "N" if lat0 >= 0 else "S"
    ew = "E" if lon0 >= 0 else "W"
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat0):02d}_00_{ew}{abs(lon0):03d}_00_DEM.tif"


tiles, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for name, kind, clat, clon in WINDOWS:
    mlat_c, mlon_c = m_per_deg(clat)
    dlat, dlon = 1.2 * HALF / mlat_c, 1.2 * HALF / mlon_c
    s = np.zeros(NF * NF)
    c = np.zeros(NF * NF)
    seen = np.zeros(NF * NF, dtype=bool)
    used, sea = [], []
    for lat0 in range(int(np.floor(clat - dlat)), int(np.floor(clat + dlat)) + 1):
        for lon0 in range(int(np.floor(clon - dlon)), int(np.floor(clon + dlon)) + 1):
            f = tiles / tile_name(lat0, lon0)
            if not f.exists():
                sea.append(f.name[22:36])
                continue
            used.append(f.name[22:36])
            z = tifffile.imread(f).astype(np.float64)
            h, w = z.shape
            lat = lat0 + 1 - np.arange(h) / h
            lon = lon0 + np.arange(w) / w
            mlat, mlon = m_per_deg(lat)
            y = (lat - clat) * mlat
            x = (lon[None, :] - clon) * mlon[:, None]
            yy = np.broadcast_to(y[:, None], (h, w))
            ix = np.floor((x + HALF) / CELL).astype(np.int64)
            iy = np.floor((yy + HALF) / CELL).astype(np.int64)
            ok = (ix >= 0) & (ix < NF) & (iy >= 0) & (iy < NF)
            k = (iy * NF + ix)[ok]
            zz = z[ok]
            seen[k] = True
            land = zz > 0.0
            np.add.at(s, k[land], zz[land])
            np.add.at(c, k[land], 1.0)
    g = np.where(c > 0, s / np.maximum(c, 1), -50.0).reshape(NF, NF)
    g = g.reshape(N, NF // N, N, NF // N).mean(axis=(1, 3))
    with open(out / f"{name}_256.bin", "wb") as fh:
        fh.write(np.uint32(N).tobytes())
        fh.write(np.uint32(N).tobytes())
        fh.write(g.astype("<f4").tobytes())
    lc = g > 0
    print(f"{name} ({kind}): tiles {used}, tiles in the margin not on disk {sea}, fine cells covered by a tile {int(seen.sum())}/{NF * NF}, land {100 * lc.mean():.1f} %, max {g.max():.0f} m")

# Corsica west: a crop of corse_1024 (195.3125 m cells; the 200 km frame centred at 42.18 N 9.045 E)
src = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("data/corsica/corse_1024.bin")
raw = np.fromfile(src, dtype="<u4", count=2)
cg = np.fromfile(src, dtype="<f4", offset=8).reshape(int(raw[1]), int(raw[0]))
mlat_c, mlon_c = m_per_deg(42.20)
ci = int(round(((8.65 - 9.045) * mlon_c + 100_000.0) / 195.3125)) - N // 2
cj = int(round(((42.20 - 42.18) * mlat_c + 100_000.0) / 195.3125)) - N // 2
g = cg[cj:cj + N, ci:ci + N]
with open(out / "corsica_west_256.bin", "wb") as fh:
    fh.write(np.uint32(N).tobytes())
    fh.write(np.uint32(N).tobytes())
    fh.write(g.astype("<f4").tobytes())
print(f"corsica_west (rocky): crop of {src.name} at columns {ci}..{ci + N}, rows {cj}..{cj + N}, land {100 * (g > 0).mean():.1f} %, max {g.max():.0f} m")
