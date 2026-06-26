"""Assemble a relief-comparison grid from per-seed images.

Usage:
    python image_grid_compare.py [<directory>] [-o output.png]

The directory is scanned for files matching `seed<NNNNN>_<col>.<ext>` where
<col> is one of altitude, boundaries, sthickness. A single grid is produced:
rows are seeds (the discriminant), columns are altitude / boundaries /
sthickness. Files that don't match the pattern are ignored.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp"}

COLUMNS = ["altitude", "boundaries", "sthickness"]

# pattern: seed<digits>_<col>
_PATTERN = re.compile(r"^seed(?P<seed>\d+)_(?P<col>altitude|boundaries|sthickness)$")

LABEL_W = 160   # left gutter for row (seed) labels
HEADER_H = 48   # top strip for column labels


def collect(directory: Path) -> dict[str, dict[str, Path]]:
    """Return {seed: {col: path}}."""
    out: dict[str, dict[str, Path]] = {}
    for path in directory.iterdir():
        if not path.is_file() or path.suffix.lower() not in IMAGE_EXTS:
            continue
        m = _PATTERN.match(path.stem)
        if not m:
            continue
        out.setdefault(m.group("seed"), {})[m.group("col")] = path
    if not out:
        raise SystemExit("error: no seed<NNNNN>_<altitude|boundaries|sthickness> images found")
    return out


def _font(size: int) -> ImageFont.ImageFont:
    for name in ("DejaVuSans.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def build_grid(seeds: dict[str, dict[str, Path]]) -> Image.Image:
    seed_keys = sorted(seeds)

    # Determine cell size from the available tiles.
    cell_w = cell_h = 0
    for cols in seeds.values():
        for p in cols.values():
            with Image.open(p) as im:
                cell_w = max(cell_w, im.width)
                cell_h = max(cell_h, im.height)
    if cell_w == 0 or cell_h == 0:
        raise SystemExit("error: no usable images found")

    grid_w = LABEL_W + cell_w * len(COLUMNS)
    grid_h = HEADER_H + cell_h * len(seed_keys)
    grid = Image.new("RGBA", (grid_w, grid_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(grid)

    header_font = _font(26)
    label_font = _font(22)

    # Column headers.
    for c, col in enumerate(COLUMNS):
        x = LABEL_W + c * cell_w + cell_w // 2
        draw.text((x, HEADER_H // 2), col, fill=(0, 0, 0, 255),
                  font=header_font, anchor="mm")

    # Rows.
    for r, seed in enumerate(seed_keys):
        y = HEADER_H + r * cell_h
        draw.text((10, y + cell_h // 2), f"seed{seed}", fill=(0, 0, 0, 255),
                  font=label_font, anchor="lm")
        for c, col in enumerate(COLUMNS):
            p = seeds[seed].get(col)
            if p is None:
                continue
            with Image.open(p) as im:
                tile = im.convert("RGBA")
            x = LABEL_W + c * cell_w + (cell_w - tile.width) // 2
            grid.paste(tile, (x, y + (cell_h - tile.height) // 2))

    return grid


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        default=Path("docs/reports/c1_continental_buoyancy/closure_morphology/relief_compare"),
        help="directory containing seed<NNNNN>_<col>.* images",
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
