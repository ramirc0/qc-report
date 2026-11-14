"""Command line: render a run's report from its metrics table and spec."""

import argparse
import os
import sys
from pathlib import Path

import polars as pl

from .render import render_report
from .spec import Spec


def build_parser():
    """Build the qc-report command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc-report flags.
    """
    p = argparse.ArgumentParser(
        prog="qc-report", description="Render a run's QC report."
    )
    p.add_argument("--spec", required=True, help="YAML spec.")
    p.add_argument(
        "--metrics",
        required=True,
        help="TSV with sample_id and one column per spec metric.",
    )
    p.add_argument(
        "--run-dir", required=True, help="Run directory. Spec links are relative to it."
    )
    p.add_argument("--run-id", required=True, help="Run name for the header.")
    p.add_argument("--version", default="", help="Pipeline version for the header.")
    p.add_argument("--output", required=True, help="HTML to write.")
    return p


def main(argv=None):
    """Render the report.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Raises
    ------
    SystemExit
        If the metrics table lacks a spec metric.
    """
    args = build_parser().parse_args(argv)
    spec = Spec.load(args.spec)
    metrics = pl.read_csv(args.metrics, separator="\t")
    ids = [metric.id for metric in spec.metrics]
    if missing := [i for i in ids if i not in metrics.columns]:
        raise SystemExit(f"{args.metrics} lacks spec metrics: {', '.join(missing)}")
    if extra := sorted(set(metrics.columns) - {"sample_id", *ids}):
        print(
            f"qc-report: not in the spec, so not shown: {', '.join(extra)}",
            file=sys.stderr,
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    run_href = Path(os.path.relpath(args.run_dir, output.parent)).as_posix()
    output.write_text(render_report(metrics, spec, args.run_id, args.version, run_href))
