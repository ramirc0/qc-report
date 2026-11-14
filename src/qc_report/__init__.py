"""One HTML QC report per run, graded by ENCODE's audits."""

from .audits import run_audits
from .cohort import find_outliers, percentiles
from .figures import (
    LEVEL_COLORS,
    gc_bias_figure,
    idr_rank_figure,
    level_histogram,
    read_quality_figure,
)
from .readers import (
    count_records,
    read_cutadapt,
    read_fastp,
    read_fingerprint,
    read_flagstat,
    read_picard_metrics,
)
from .render import format_value, render_report, write_sample_panels
from .spec import Audit, Group, Level, Link, Metric, Spec, SpecError, Standards

__all__ = [
    "LEVEL_COLORS",
    "Audit",
    "Group",
    "Level",
    "Link",
    "Metric",
    "Spec",
    "SpecError",
    "Standards",
    "count_records",
    "find_outliers",
    "format_value",
    "gc_bias_figure",
    "idr_rank_figure",
    "level_histogram",
    "percentiles",
    "read_cutadapt",
    "read_fastp",
    "read_fingerprint",
    "read_flagstat",
    "read_picard_metrics",
    "read_quality_figure",
    "render_report",
    "run_audits",
    "write_sample_panels",
]
