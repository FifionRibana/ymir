"""Assemble a climate-verdict grid from per-seed images, with labelled legends.

Usage:
    python image_grid_climate.py [<directory>] [-o output.png]

The directory is scanned for files matching `seed<NNNNN>_<...>.<ext>`. Each seed
contributes one row; columns are relief / precip / biomes / temperature /
transect (matched by keyword, the variable suffixes like `_wind3577m_HAUTE` or
`_row550` are ignored for classification). Three labelled legends are drawn at
the bottom:

  - precip bands  (legend_precip_bands.png): 5 mm/yr slices, beige -> green ->
    blue = dry -> wet. NOT a biome palette.
  - biomes        (legend_biomes.png): 6-colour categorical palette.
  - temp bands    (legend_temp_bands.png): 4 slices at the Whittaker thresholds
    -5 / +5 / +20 C, cold -> hot. Identical scale across all seeds.

The legend colours/labels are drawn from the tables below (the source legend
PNGs are bare colour strips with no text).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp"}

# Columns: (display name, keyword used to classify the file's suffix).
COLUMNS = [
    ("relief", "relief"),
    ("precip", "precip"),
    ("biomes", "biomes"),
    ("temperature", "temp"),
    ("transect", "transect"),
]

# Precip bands: 5 mm/yr slices, dry -> wet (beige -> green -> blue).
PRECIP_LEGEND = [
    ((225, 200, 140), "desert  < 250 mm/an"),
    ((200, 195, 110), "steppe  250-500"),
    ((150, 180, 90), "tempere-sec  500-800"),
    ((80, 150, 200), "oceanique  800-1500"),
    ((30, 90, 200), "tres humide  > 1500"),
]

# Biomes: 6-colour categorical palette.
BIOME_LEGEND = [
    ((200, 195, 110), "steppe (temperate grassland)"),
    ((80, 160, 80), "foret temperee"),
    ((40, 110, 70), "foret temperee humide"),
    ((70, 110, 90), "boreal / taiga"),
    ((200, 205, 215), "toundra"),
    ((225, 200, 140), "desert"),
]

# Temperature: 4 bands at the Whittaker thresholds -5 / +5 / +20 C (cold -> hot).
TEMP_LEGEND = [
    ((225, 235, 248), "polaire  < -5 C"),
    ((90, 140, 205), "boreal  -5 a +5"),
    ((110, 190, 110), "tempere  +5 a +20"),
    ((225, 120, 70), "chaud  > +20"),
]

# pattern: seed<digits>_<rest>
_PATTERN = re.compile(r"^seed(?P<seed>\d+)_(?P<rest>.+)$")

LABEL_W = 160    # left gutter for row (seed) labels
HEADER_H = 48    # top strip for column labels
LEGEND_H = 360   # bottom strip for the two legends
PAD = 12
ROW_H = 480      # tiles are scaled to this height (aspect preserved)


def collect(directory: Path) -> dict[str, dict[str, Path]]:
    """Return {seed: {column-name: path}}."""
    out: dict[str, dict[str, Path]] = {}
    for path in directory.iterdir():
        if not path.is_file() or path.suffix.lower() not in IMAGE_EXTS:
            continue
        if path.stem.startswith("legend"):
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
        raise SystemExit("error: no seed<NNNNN>_... images found")
    return out


def _font(size: int) -> ImageFont.ImageFont:
    for name in ("DejaVuSans.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def draw_legend(draw: ImageDraw.ImageDraw, x: int, y: int, title: str,
                entries: list[tuple[tuple[int, int, int], str]],
                title_font: ImageFont.ImageFont, font: ImageFont.ImageFont) -> None:
    draw.text((x, y), title, fill=(0, 0, 0, 255), font=title_font, anchor="lm")
    sw = 34   # swatch width
    sh = 24   # swatch height
    row_h = 34
    yy = y + 28
    for color, label in entries:
        draw.rectangle([x, yy, x + sw, yy + sh], fill=color + (255,),
                       outline=(60, 60, 60, 255))
        draw.text((x + sw + 12, yy + sh // 2), label, fill=(0, 0, 0, 255),
                  font=font, anchor="lm")
        yy += row_h


def _load_scaled(path: Path) -> Image.Image:
    """Open a tile and scale it to ROW_H height, preserving aspect ratio."""
    im = Image.open(path).convert("RGBA")
    if im.height != ROW_H:
        w = max(1, round(im.width * ROW_H / im.height))
        im = im.resize((w, ROW_H), Image.LANCZOS)
    return im


def build_grid(seeds: dict[str, dict[str, Path]]) -> Image.Image:
    seed_keys = sorted(seeds)
    col_names = [name for name, _ in COLUMNS]

    # Scale every tile to a common row height; column width = widest tile there.
    row_h = ROW_H
    col_w = {name: 0 for name in col_names}
    for cols in seeds.values():
        for name, p in cols.items():
            col_w[name] = max(col_w[name], _load_scaled(p).width)
    if all(w == 0 for w in col_w.values()):
        raise SystemExit("error: no usable images found")

    col_x = {}
    x_acc = LABEL_W
    for name in col_names:
        col_x[name] = x_acc
        x_acc += col_w[name]

    grid_w = x_acc
    grid_h = HEADER_H + row_h * len(seed_keys) + LEGEND_H
    grid = Image.new("RGBA", (grid_w, grid_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(grid)

    header_font = _font(28)
    label_font = _font(22)
    legend_title_font = _font(24)
    legend_font = _font(20)

    # Column headers.
    for name in col_names:
        x = col_x[name] + col_w[name] // 2
        draw.text((x, HEADER_H // 2), name, fill=(0, 0, 0, 255),
                  font=header_font, anchor="mm")

    # Rows (seeds).
    for r, seed in enumerate(seed_keys):
        y = HEADER_H + r * row_h
        draw.text((10, y + row_h // 2), f"seed{seed}", fill=(0, 0, 0, 255),
                  font=label_font, anchor="lm")
        for name in col_names:
            p = seeds[seed].get(name)
            if p is None:
                continue
            tile = _load_scaled(p)
            x = col_x[name] + (col_w[name] - tile.width) // 2
            grid.paste(tile, (x, y + (row_h - tile.height) // 2))

    # Legends at the bottom: three columns (precip / biomes / temperature).
    legend_y = HEADER_H + row_h * len(seed_keys) + PAD
    draw.line([(0, legend_y - PAD // 2), (grid_w, legend_y - PAD // 2)],
              fill=(180, 180, 180, 255), width=2)
    legends = [
        ("bandes de precipitation (mm/an, sec -> humide)", PRECIP_LEGEND),
        ("biomes (palette categorielle)", BIOME_LEGEND),
        ("temperature (C, seuils Whittaker -5/+5/+20)", TEMP_LEGEND),
    ]
    usable_w = grid_w - LABEL_W
    for i, (title, entries) in enumerate(legends):
        x = LABEL_W + PAD + i * usable_w // len(legends)
        draw_legend(draw, x, legend_y + 14, title, entries,
                    legend_title_font, legend_font)

    return grid


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        default=Path("docs/reports/c1_continental_buoyancy/closure_morphology/climate_verdict"),
        help="directory containing seed<NNNNN>_... images",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="output image path (default: <directory>/grid.png)",
    )
    args = parser.parse_args()

    if not args.directory.is_dir():
        raise SystemExit(f"error: not a directory: {args.directory}")

    output = args.output or args.directory / "grid.png"
    seeds = collect(args.directory)
    grid = build_grid(seeds)
    output.parent.mkdir(parents=True, exist_ok=True)
    grid.save(output)
    print(f"wrote {output} ({grid.width}x{grid.height}, {len(seeds)} seeds)")


if __name__ == "__main__":
    main()
