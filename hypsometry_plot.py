"""Plot the C1 hypsometric diagnostic (#165 amont) vs Earth's land hypsography.

Reads hypsometry_bands.csv (written by the `probe_hypsometry` Rust test) and
produces two PNGs in the same directory:
  - hypsometry_bands.png      : band-fraction bars, all/coast/interior vs Earth.
  - hypsometry_cumulative.png : cumulative "fraction of land below X metres".
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = Path("docs/reports/c1_continental_buoyancy/closure_morphology/hypsometry")
BANDS = ["0-250", "250-500", "500-1000", "1000-2000", "2000-3000", "3000+"]
UPPER = [250, 500, 1000, 2000, 3000, 4000]  # band upper edges (m); last capped


def load() -> dict[tuple[str, str], list[float]]:
    """Return {(seed, region): [frac per band]} for AGGREGATE + EARTH rows."""
    out: dict[tuple[str, str], dict[str, float]] = {}
    with open(DIR / "hypsometry_bands.csv", newline="") as f:
        for row in csv.DictReader(f):
            key = (row["seed"], row["region"])
            out.setdefault(key, {})[row["band_m"]] = float(row["fraction"])
    return {k: [v[b] for b in BANDS] for k, v in out.items()}


def cumulative(fracs: list[float]) -> list[float]:
    acc, out = 0.0, []
    for v in fracs:
        acc += v
        out.append(acc)
    return out


def main() -> None:
    data = load()
    agg_all = data[("AGGREGATE", "all")]
    agg_coast = data[("AGGREGATE", "coast")]
    agg_int = data[("AGGREGATE", "interior")]
    earth = data[("EARTH", "land")]

    # --- 1. band-fraction bars ---
    fig, ax = plt.subplots(figsize=(11, 6))
    x = range(len(BANDS))
    width = 0.2
    ax.bar([i - 1.5 * width for i in x], earth, width, label="Earth (land)", color="#888888")
    ax.bar([i - 0.5 * width for i in x], agg_all, width, label="C1 all land", color="#4c72b0")
    ax.bar([i + 0.5 * width for i in x], agg_coast, width, label="C1 coast (<=150km)", color="#55a868")
    ax.bar([i + 1.5 * width for i in x], agg_int, width, label="C1 interior (>150km)", color="#c44e52")
    ax.set_xticks(list(x))
    ax.set_xticklabels(BANDS)
    ax.set_xlabel("altitude band (m)")
    ax.set_ylabel("fraction of land area")
    ax.set_title("#165 hypsometry — land altitude distribution vs Earth (6 seeds aggregated)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(DIR / "hypsometry_bands.png", dpi=130)
    plt.close(fig)

    # --- 2. cumulative "fraction below X m" ---
    fig, ax = plt.subplots(figsize=(11, 6))
    for label, fr, color in [
        ("Earth (land)", earth, "#888888"),
        ("C1 all land", agg_all, "#4c72b0"),
        ("C1 coast (<=150km)", agg_coast, "#55a868"),
        ("C1 interior (>150km)", agg_int, "#c44e52"),
    ]:
        ax.plot(UPPER, cumulative(fr), marker="o", label=label, color=color)
    ax.axvline(900, ls="--", color="orange", lw=1.5,
               label="~900 m (temperate-forest ceiling @45°)")
    ax.set_xlabel("altitude (m)")
    ax.set_ylabel("fraction of land below altitude")
    ax.set_title("#165 hypsometry — cumulative curve (fraction of land below X m)")
    ax.set_ylim(0, 1.02)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(DIR / "hypsometry_cumulative.png", dpi=130)
    plt.close(fig)

    print(f"wrote {DIR / 'hypsometry_bands.png'} and {DIR / 'hypsometry_cumulative.png'}")


if __name__ == "__main__":
    main()
