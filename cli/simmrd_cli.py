#!/usr/bin/env python3
"""simmrd-cli: run simmrd simulations from a YAML parameter file.

Usage (inside the pixi environment):
    python simmrd_cli.py --params params/example_summary.yaml --output results.rds
    pixi run simulate  --params params/example_summary.yaml --output results.rds --seed 42
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="simmrd",
        description="Run simmrd simulations defined by a YAML parameter file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  simmrd --params params/example_summary.yaml --output out.rds\n"
            "  simmrd --params params/example_summary.yaml --output out.rds -n 500\n"
            "  simmrd --params params/example_summary.yaml --output out.rds -n 500 --seed 42\n"
            "\n"
            "See params/ for annotated example YAML files.\n"
            "Run 'pixi run setup' once to install simmrd from source before first use.\n"
        ),
    )
    parser.add_argument(
        "--params",
        required=True,
        metavar="FILE",
        help="YAML file containing simmrd parameters (see params/ for examples).",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Destination path for the output dataset (format TBD; currently RDS).",
    )
    parser.add_argument(
        "-n", "--iterations",
        type=int,
        default=1,
        metavar="INT",
        help="Number of simulation iterations (default: 1). "
             "Results are stored as a list and saved to --output.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="Master integer seed. When set, all iterations are fully "
             "reproducible. Omit for a fresh random sequence each run.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.iterations < 1:
        sys.exit("Error: --iterations must be >= 1")

    runner = HERE / "runner.R"
    if not runner.exists():
        sys.exit(f"Error: runner.R not found at {runner}")

    cmd: list[str] = [
        "Rscript", str(runner),
        "--params", args.params,
        "--output", args.output,
        "--iterations", str(args.iterations),
    ]
    if args.seed is not None:
        cmd += ["--seed", str(args.seed)]

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
