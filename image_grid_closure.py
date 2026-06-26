"""Assemble per-seed closure-morphology grids from mode/seed images.

Usage:
    python image_grid_closure.py <directory> [-o output_dir]

The directory is scanned for files matching `<mode>_seed<NNNNN>_altitude.<ext>`
and `<mode>_seed<NNNNN>_s.<ext>`. One grid is produced per seed (the seed is
the discriminant): columns are `altitude` then `s`, rows are the closure modes
in a fixed order (full, no_davis_suppe, no_erosion, no_stein_stein,
no_subduction_accretion, no_track_d). Files that don't match the
`<mode>_seed<NNNNN>_<col>` pattern (e.g. trajectory frames) are ignored.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp"}

# Row order for the grid (modes). Modes found on disk but not listed here are
# appended after these, alphabetically.
MODE_ORDER = [
    "full",
    "no_davis_suppe",
    "no_erosion",
    "no_stein_stein",
    "no_subduction_accretion",
    "no_track_d",
]

COLUMNS = ["altitude", "s"]

# pattern: <mode>_seed<digits>_<col>
_PATTERN = re.compile(r"^(?P<mode>.+)_seed(?P<seed>\d+)_(?P<col>altitude|s)$")

LABEL_W = 220   # left gutter for row (mode) labels
HEADER_H = 48   # top strip for column labels


def collect(directory: Path) -> dict[str, dict[str, dict[str, Path]]]:
    """Return {seed: {mode: {col: path}}}."""
    out: dict[str, dict[str, dict[str, Path]]] = {}
    for path in directory.iterdir():
        if not path.is_file() or path.suffix.lower() not in IMAGE_EXTS:
            continue
        m = _PATTERN.match(path.stem)
        if not m:
            continue
        seed = m.group("seed")
        mode = m.group("mode")
        col = m.group("col")
        out.setdefault(seed, {}).setdefault(mode, {})[col] = path
    if not out:
        raise SystemExit("error: no <mode>_seed<NNNNN>_<altitude|s> images found")
    return out


def order_modes(modes: set[str]) -> list[str]:
    known = [m for m in MODE_ORDER if m in modes]
    extra = sorted(modes - set(MODE_ORDER))
    return known + extra


def _font(size: int) -> ImageFont.ImageFont:
    for name in ("DejaVuSans.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def build_seed_grid(seed: str, modes_map: dict[str, dict[str, Path]]) -> Image.Image:
    modes = order_modes(set(modes_map))

    # Determine cell size from the available tiles.
    cell_w = cell_h = 0
    for mode in modes:
        for col in COLUMNS:
            p = modes_map[mode].get(col)
            if p is None:
                continue
            with Image.open(p) as im:
                cell_w = max(cell_w, im.width)
                cell_h = max(cell_h, im.height)
    if cell_w == 0 or cell_h == 0:
        raise SystemExit(f"error: no usable images for seed {seed}")

    grid_w = LABEL_W + cell_w * len(COLUMNS)
    grid_h = HEADER_H + cell_h * len(modes)
    grid = Image.new("RGBA", (grid_w, grid_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(grid)

    header_font = _font(26)
    label_font = _font(20)

    # Column headers.
    for c, col in enumerate(COLUMNS):
        x = LABEL_W + c * cell_w + cell_w // 2
        draw.text((x, HEADER_H // 2), col, fill=(0, 0, 0, 255),
                  font=header_font, anchor="mm")

    # Rows.
    for r, mode in enumerate(modes):
        y = HEADER_H + r * cell_h
        draw.text((10, y + cell_h // 2), mode, fill=(0, 0, 0, 255),
                  font=label_font, anchor="lm")
        for c, col in enumerate(COLUMNS):
            p = modes_map[mode].get(col)
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
        default=Path("docs/reports/c1_continental_buoyancy/closure_morphology"),
        help="directory containing <mode>_seed<NNNNN>_<altitude|s>.* images",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="output directory (default: <directory>); writes grid_seed<NNNNN>.png",
    )
    args = parser.parse_args()

    if not args.directory.is_dir():
        raise SystemExit(f"error: not a directory: {args.directory}")

    out_dir = args.output or args.directory
    out_dir.mkdir(parents=True, exist_ok=True)

    seeds = collect(args.directory)
    for seed in sorted(seeds):
        grid = build_seed_grid(seed, seeds[seed])
        output = out_dir / f"grid_seed{seed}.png"
        grid.save(output)
        n_rows = len(order_modes(set(seeds[seed])))
        print(f"wrote {output} ({grid.width}x{grid.height}, {n_rows} modes)")


if __name__ == "__main__":
    main()
