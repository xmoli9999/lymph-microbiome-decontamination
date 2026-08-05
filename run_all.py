#!/usr/bin/env python3
"""Run the full independent analysis, acceptance checks, source-data export and Figure 2-5 rendering."""

from __future__ import annotations
import argparse
from pathlib import Path
import independent_rewrite_analysis as analysis
import render_figures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("otu_workbook", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    input_path = args.otu_workbook.resolve()
    output_dir = args.output_dir.resolve()
    analysis.run(input_path, output_dir)
    render_figures.render_all(input_path, output_dir)
    print(f"\nComplete analysis and figures written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
