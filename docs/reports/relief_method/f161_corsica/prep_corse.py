"""ADR Finding 160->161 (F161-C) -- the Corsica reference grids from Copernicus DEM GLO-30.

Usage: python -I prep_corse.py <tile dir> <out dir>

* reads the five 1-arcsecond tiles (PixelIsPoint: pixel (i, j) centre at lon0 + i/3600, lat0 + 1 - j/3600);
* projects every pixel centre to a local metric frame centred on Corsica (lat0 42.18, lon0 9.045; metres per degree
  of latitude and longitude from the WGS84 ellipsoid at the pixel's own latitude);
* averages the pixels into a 4096 x 4096 base grid of 48.828125 m cells (= 400 km / 8192, the cascade's finest cell),
  200 km on a side, rows SOUTH-FIRST (Ymir's convention);
* the sea: a base cell with no land pixel (DEM <= 0) is set to -50 m; land keeps its mean over its land pixels;
* block means to 2048 / 1024 / 512 / 256 / 128 cells (97.7 / 195.3 / 390.6 / 781.3 / 1562.5 m);
* writes corse_<n>.bin (u32 w, u32 h, f32 LE, south-first) for n in 4096 ... 128.
"""
import math
import sys
from pathlib import Path

import numpy as np
import tifffile

tiles, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
LAT0, LON0 = 42.18, 9.045
N, CELL = 4096, 400_000.0 / 8192
HALF = N * CELL / 2


def m_per_deg(lat):
    p = np.radians(lat)
    mlat = 111132.954 - 559.822 * np.cos(2 * p) + 1.175 * np.cos(4 * p)
    mlon = 111412.84 * np.cos(p) - 93.5 * np.cos(3 * p)
    return mlat, mlon


s = np.zeros(N * N)
c = np.zeros(N * N)
any_px = np.zeros(N * N, dtype=bool)
for f in sorted(tiles.glob("Copernicus_DSM_COG_10_*_DEM.tif")):
    name = f.name.split("_")
    lat0 = int(name[4][1:])
    lon0 = int(name[6][1:])
    z = tifffile.imread(f).astype(np.float64)
    h, w = z.shape
    lat = lat0 + 1 - np.arange(h) / h  # row j
    lon = lon0 + np.arange(w) / w  # column i
    mlat, mlon = m_per_deg(lat)
    y = (lat - LAT0) * mlat  # per row
    x = (lon[None, :] - LON0) * mlon[:, None]  # per row and column
    yy = np.broadcast_to(y[:, None], (h, w))
    ix = np.floor((x + HALF) / CELL).astype(np.int64)
    iy = np.floor((yy + HALF) / CELL).astype(np.int64)
    ok = (ix >= 0) & (ix < N) & (iy >= 0) & (iy < N)
    k = (iy * N + ix)[ok]
    zz = z[ok]
    any_px[k] = True
    land = zz > 0.0
    np.add.at(s, k[land], zz[land])
    np.add.at(c, k[land], 1.0)
    print(f.name, z.shape, "pixels in the grid:", int(ok.sum()), "land:", int(land.sum()), "max", float(z.max()))
base = np.where(c > 0, s / np.maximum(c, 1), -50.0).reshape(N, N)
# keep Corsica alone: the largest 4-connected land mass; the other islands in the frame (Capraia, La Maddalena) -> sea
land = base > 0
lab = np.zeros(N * N, dtype=np.int32)
flat = land.ravel()
sizes = {}
cur = 0
for s0 in np.flatnonzero(flat):
    if lab[s0]:
        continue
    cur += 1
    stack = [s0]
    lab[s0] = cur
    n_ = 0
    while stack:
        k = stack.pop()
        n_ += 1
        x, y = k % N, k // N
        for j in ((k - 1) if x > 0 else -1, (k + 1) if x + 1 < N else -1, (k - N) if y > 0 else -1, (k + N) if y + 1 < N else -1):
            if j >= 0 and flat[j] and not lab[j]:
                lab[j] = cur
                stack.append(j)
    sizes[cur] = n_
main = max(sizes, key=sizes.get)
print("land masses:", len(sizes), "Corsica", sizes[main], "cells; dropped", sum(v for k_, v in sizes.items() if k_ != main), "cells")
base = np.where(lab.reshape(N, N) == main, base, -50.0)
print("base: land cells", int((c > 0).sum()), "max", float(base.max()), "cells without any DEM pixel", int((~any_px).sum()))


def save(g, n):
    with open(out / f"corse_{n}.bin", "wb") as fh:
        fh.write(np.uint32(n).tobytes())
        fh.write(np.uint32(n).tobytes())
        fh.write(g.astype("<f4").tobytes())


g = base
n = N
save(g, n)
while n > 128:
    g = g.reshape(n // 2, 2, n // 2, 2).mean(axis=(1, 3))
    n //= 2
    save(g, n)
    land = g > 0
    print(f"corse_{n}: cell {200_000 / n:.1f} m, land cells {int(land.sum())}, max {float(g.max()):.0f} m, mean land {float(g[land].mean()):.0f} m")
