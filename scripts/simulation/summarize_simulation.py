from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = pd.read_csv(args.input)
    summary = data[["events", "valid_lower", "valid_upper"]].agg(["count", "mean", "std"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output)


if __name__ == "__main__":
    main()
