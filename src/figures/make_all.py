"""
Render every report figure into outputs/report_figures/.

Assumes outputs/report_figures/_stats.json already exists; build it first with
    python src/figures/compute_stats.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fig_architecture
import fig_method
import fig_reach
import fig_validation
from common import OUT, STATS


def main():
    if not STATS.exists():
        sys.exit(f"missing {STATS} — run compute_stats.py first")
    for mod in (fig_architecture, fig_reach, fig_method, fig_validation):
        mod.main()
    pngs = sorted(OUT.glob("*.png"))
    print(f"\n{len(pngs)} figures in {OUT}")
    for p in pngs:
        print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
