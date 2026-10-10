"""F165-B1: ETOPO 2022 15" surface tiles -> Europe (3 000 km) and an alpine window (400 km), Lambert azimuthal
equal-area (spherical, R = 6 371 km), 1.5625 km cells, each the mean of 4x4 bilinear samples. Rows south-first,
u32 w, u32 h, f32 LE (the data/corsica convention). Also a region-id grid for the qualitative control.
Usage: python -I prep_europe.py <tiles dir> <out dir>"""
import sys, os, glob
import numpy as np
import tifffile

tiles_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
R = 6371.0
CELL = 1.5625
SUB = 4

# mosaic lat 30..75, lon -15..45 at 15" (240 px per degree)
PPD = 240
LAT_TOP, LON_LEFT = 75.0, -15.0
H, W = 45 * PPD, 60 * PPD
mosaic = np.full((H, W), np.nan, dtype=np.float32)
for p in sorted(glob.glob(os.path.join(tiles_dir, 'ETOPO_2022_v1_15s_*_surface.tif'))):
    with tifffile.TiffFile(p) as t:
        pg = t.pages[0]
        tie = pg.tags['ModelTiepointTag'].value
        a = pg.asarray().astype(np.float32)
    lon0, lat0 = tie[3], tie[4]
    r0 = int(round((LAT_TOP - lat0) * PPD))
    c0 = int(round((lon0 - LON_LEFT) * PPD))
    a[a <= -99998] = np.nan
    mosaic[r0:r0 + a.shape[0], c0:c0 + a.shape[1]] = a
    print('tile', os.path.basename(p), 'at', r0, c0, a.shape)

def sample(lat, lon):
    fr = (LAT_TOP - lat) * PPD - 0.5
    fc = (lon - LON_LEFT) * PPD - 0.5
    r = np.floor(fr).astype(np.int64); c = np.floor(fc).astype(np.int64)
    tr = fr - r; tc = fc - c
    ok = (r >= 0) & (c >= 0) & (r + 1 < H) & (c + 1 < W)
    r = np.clip(r, 0, H - 2); c = np.clip(c, 0, W - 2)
    v = (mosaic[r, c] * (1 - tr) * (1 - tc) + mosaic[r, c + 1] * (1 - tr) * tc
         + mosaic[r + 1, c] * tr * (1 - tc) + mosaic[r + 1, c + 1] * tr * tc)
    v = np.where(ok, v, np.nan)
    return v

def laea_inv(x, y, lat1, lon1):
    p1, l0 = np.radians(lat1), np.radians(lon1)
    rho = np.hypot(x, y)
    c = 2 * np.arcsin(np.clip(rho / (2 * R), -1, 1))
    with np.errstate(invalid='ignore', divide='ignore'):
        lat = np.arcsin(np.cos(c) * np.sin(p1) + np.where(rho > 0, y * np.sin(c) * np.cos(p1) / rho, 0))
        lon = l0 + np.arctan2(x * np.sin(c), rho * np.cos(p1) * np.cos(c) - y * np.sin(p1) * np.sin(c))
    return np.degrees(lat), np.degrees(lon)

def window(lat1, lon1, n):
    half = n * CELL / 2
    out = np.zeros((n, n), dtype=np.float64)
    cnt = np.zeros((n, n), dtype=np.float64)
    missing = np.zeros((n, n), dtype=np.int64)
    lat_c = np.zeros((n, n)); lon_c = np.zeros((n, n))
    for row in range(0, n, 64):
        rows = np.arange(row, min(row + 64, n))
        for si in range(SUB):
            for sj in range(SUB):
                xs = -half + (np.arange(n) + (sj + 0.5) / SUB) * CELL
                ys = -half + (rows + (si + 0.5) / SUB) * CELL   # south-first: row 0 is the south
                X, Y = np.meshgrid(xs, ys)
                la, lo = laea_inv(X, Y, lat1, lon1)
                v = sample(la, lo)
                good = np.isfinite(v)
                out[rows] += np.where(good, v, 0)
                cnt[rows] += good
                missing[rows] += ~good
        xs = -half + (np.arange(n) + 0.5) * CELL
        ys = -half + (rows + 0.5) * CELL
        X, Y = np.meshgrid(xs, ys)
        la, lo = laea_inv(X, Y, lat1, lon1)
        lat_c[rows] = la; lon_c[rows] = lo
    z = np.where(cnt > 0, out / np.maximum(cnt, 1), -3000.0).astype(np.float32)
    return z, lat_c, lon_c, int((cnt == 0).sum())

def save(path, z):
    n = z.shape[0]
    with open(path, 'wb') as f:
        f.write(np.uint32(z.shape[1]).tobytes()); f.write(np.uint32(z.shape[0]).tobytes())
        f.write(z.astype('<f4').tobytes())

z, la, lo, miss = window(50.0, 10.0, 1920)
print('europe: no-data cells set to -3000 m:', miss, 'lat', la.min(), la.max(), 'lon', lo.min(), lo.max())
save(os.path.join(out_dir, 'europe_1920.bin'), z)
reg = np.zeros_like(z)
boxes = {1: (47.8, 49.3, 1.0, 3.5), 2: (44.7, 46.0, 2.3, 4.0), 3: (39.0, 41.5, -5.5, -2.5), 4: (41.4, 43.0, 8.5, 9.6)}
for k, (a, b, c, d) in boxes.items():
    reg[(la >= a) & (la <= b) & (lo >= c) & (lo <= d)] = k
save(os.path.join(out_dir, 'europe_1920_regions.bin'), reg)
print('regions cells', {k: int((reg == k).sum()) for k in boxes})
za, la2, lo2, miss2 = window(46.3, 9.5, 256)
print('alps: missing', miss2, 'lat', la2.min(), la2.max(), 'lon', lo2.min(), lo2.max(), 'max z', za.max())
save(os.path.join(out_dir, 'alps_256.bin'), za)
print('europe max z', float(z.max()), 'land share', float((z > 0).mean()))
