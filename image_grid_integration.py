"""#viz — assemble the INTEGRATION grid: the whole current chain on one shared
terrain per seed (relief+bathymetry / drainage / precip / temperature / biomes /
transect), 6 seeds, common scale. The global visual-validation instrument before
the Living Landz export.

Tiles come from `probe_integration_grid`
(`docs/reports/c1_continental_buoyancy/closure_morphology/integration_grid/`).
This only assembles + labels; the Rust probe bakes the common scales.

Usage:
    python image_grid_integration.py [<dir>] [-o out.png]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

TILES_DIR = Path("docs/reports/c1_continental_buoyancy/closure_morphology/integration_grid")

# (display name, filename keyword). Suffixes (_MMBANDS/_CBANDS) ignored.
COLUMNS = [
    ("relief + bathy", "relief"),
    ("drainage", "drainage"),
    ("precip", "precip"),
    ("temperature", "temp"),
    ("biomes", "biomes"),
    ("transect", "transect"),
]

RELIEF_LEGEND = [
    ((10, 25, 70), "abysse  ~ −5000 m"),
    ((140, 195, 215), "plateau  ~ 0 m"),
    ((70, 130, 70), "plaine  0-1000 m"),
    ((200, 180, 120), "collines  1000-3000 m"),
    ((230, 230, 235), "sommets  > 3000 m"),
]
DRAINAGE_LEGEND = [
    ((90, 160, 240), "rivière  small-boat"),
    ((40, 110, 230), "rivière  barge"),
    ((20, 70, 200), "rivière  ship"),
    ((30, 90, 180), "lac"),
]
PRECIP_LEGEND = [
    ((225, 200, 140), "desert  < 250 mm/an"),
    ((200, 195, 110), "steppe  250-500"),
    ((150, 180, 90), "tempéré-sec  500-800"),
    ((80, 150, 200), "océanique  800-1500"),
    ((30, 90, 200), "très humide  > 1500"),
]
TEMP_LEGEND = [
    ((225, 235, 248), "polaire  < −5 °C"),
    ((90, 140, 205), "boréal  −5 à +5"),
    ((110, 190, 110), "tempéré  +5 à +20"),
    ((225, 120, 70), "chaud  > +20"),
]
BIOME_LEGEND = [
    ((200, 195, 110), "steppe"),
    ((80, 160, 80), "forêt tempérée"),
    ((40, 110, 70), "forêt temp. humide"),
    ((70, 110, 90), "boréal / taïga"),
    ((200, 205, 215), "toundra"),
    ((20, 110, 50), "forêt tropicale"),
]

_PATTERN = re.compile(r"^seed(?P<seed>\d+)_(?P<rest>.+)$")
LABEL_W = 130
HEADER_H = 44
LEGEND_H = 250
PAD = 12
ROW_H = 430


def collect(directory: Path) -> dict[str, dict[str, Path]]:
    out: dict[str, dict[str, Path]] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() != ".png" or path.stem.startswith("legend"):
            continue
        m = _PATTERN.match(path.stem)
        if not m:
            continue
        rest = m.group("rest").lower()
        for name, keyword in COLUMNS:
            if keyword in rest:
                out.setdefault(m.group("seed"), {})[name] = path
                break
    if not out:
        raise SystemExit(f"error: no seed<NNNNN>_*.png in {directory} (run probe_integration_grid)")
    return out


def font(size: int) -> ImageFont.ImageFont:
    for name in ("DejaVuSans.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def load_scaled(path: Path) -> Image.Image:
    im = Image.open(path).convert("RGBA")
    if im.height != ROW_H:
        w = max(1, round(im.width * ROW_H / im.height))
        im = im.resize((w, ROW_H), Image.LANCZOS)
    return im


def draw_legend(draw, x, y, title, entries, tf, lf):
    draw.text((x, y), title, fill=(0, 0, 0, 255), font=tf, anchor="lm")
    sw, sh, rh = 26, 18, 26
    yy = y + 24
    for color, label in entries:
        draw.rectangle([x, yy, x + sw, yy + sh], fill=color + (255,), outline=(60, 60, 60, 255))
        draw.text((x + sw + 8, yy + sh // 2), label, fill=(0, 0, 0, 255), font=lf, anchor="lm")
        yy += rh


def build(seeds: dict[str, dict[str, Path]]) -> Image.Image:
    keys = sorted(seeds)
    names = [n for n, _ in COLUMNS]
    col_w = {n: 0 for n in names}
    for cols in seeds.values():
        for n in names:
            if n in cols:
                col_w[n] = max(col_w[n], load_scaled(cols[n]).width)
    col_x, x = {}, LABEL_W
    for n in names:
        col_x[n] = x
        x += col_w[n]
    grid_w, grid_h = x, HEADER_H + ROW_H * len(keys) + LEGEND_H
    grid = Image.new("RGBA", (grid_w, grid_h), (255, 255, 255, 255))
    d = ImageDraw.Draw(grid)
    fh, fl, ft, flg = font(24), font(20), font(22), font(17)

    for n in names:
        d.text((col_x[n] + col_w[n] // 2, HEADER_H // 2), n, fill=(0, 0, 0, 255), font=fh, anchor="mm")
    for r, seed in enumerate(keys):
        y = HEADER_H + r * ROW_H
        d.text((8, y + ROW_H // 2), f"seed{int(seed)}", fill=(0, 0, 0, 255), font=fl, anchor="lm")
        for n in names:
            p = seeds[seed].get(n)
            if p is None:
                continue
            t = load_scaled(p)
            grid.paste(t, (col_x[n] + (col_w[n] - t.width) // 2, y + (ROW_H - t.height) // 2))

    ly = HEADER_H + ROW_H * len(keys) + PAD
    d.line([(0, ly - PAD // 2), (grid_w, ly - PAD // 2)], fill=(180, 180, 180, 255), width=2)
    legends = [
        ("relief + bathymétrie", RELIEF_LEGEND),
        ("drainage", DRAINAGE_LEGEND),
        ("précipitation (mm/an)", PRECIP_LEGEND),
        ("température (°C)", TEMP_LEGEND),
        ("biomes", BIOME_LEGEND),
    ]
    usable = grid_w - LABEL_W
    for i, (title, entries) in enumerate(legends):
        draw_legend(d, LABEL_W + PAD + i * usable // len(legends), ly + 14, title, entries, ft, flg)
    return grid


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory", nargs="?", type=Path, default=TILES_DIR)
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()
    if not args.directory.is_dir():
        raise SystemExit(f"error: not a directory: {args.directory} (run probe_integration_grid first)")
    seeds = collect(args.directory)
    grid = build(seeds)
    out = args.output or args.directory / "grid_integration.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    grid.save(out)
    print(f"wrote {out} ({grid.width}x{grid.height}, {len(seeds)} seeds)")
    print("COMMON SCALE baked by probe_integration_grid (relief/bathy, precip, temp, biomes) — comparable across seeds.")


if __name__ == "__main__":
    main()
