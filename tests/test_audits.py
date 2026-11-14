import polars as pl
import pytest

from qc_report import Audit, Group, Level, Metric, Spec, run_audits
from qc_report.audits import audit_bounds, cell_levels, sample_levels


def atac_spec():
    metrics = [
        Metric("alignment_pct", "Aligned", format="percent"),
        Metric("overlap_peaks", "Overlap peaks", format="count"),
        Metric("idr_peaks", "IDR peaks", format="count"),
        Metric("nfr_peak", "NFR peak", format="flag"),
    ]
    audits = [
        Audit(
            "acceptable alignment rate",
            "WARNING",
            "alignment_pct >= 80 AND alignment_pct <= 95",
        ),
        Audit("low alignment rate", "NOT_COMPLIANT", "alignment_pct < 80"),
        Audit(
            "moderate number of reproducible peaks",
            "WARNING",
            "overlap_peaks <= 150000 AND idr_peaks <= 70000"
            " AND (overlap_peaks >= 100000 OR idr_peaks >= 50000)",
        ),
        Audit(
            "insufficient number of reproducible peaks",
            "NOT_COMPLIANT",
            "overlap_peaks < 100000 AND idr_peaks < 50000",
        ),
        Audit("no peaks in nucleosome-free regions", "NOT_COMPLIANT", "NOT nfr_peak"),
    ]
    return Spec("ATAC", [Group("All", metrics)], audits)


def levels_of(rows):
    metrics = pl.DataFrame(rows)
    status = sample_levels(run_audits(metrics, atac_spec()))
    return [status.get(s, Level.OK) for s in metrics["sample_id"]]


@pytest.mark.parametrize(
    ("aligned", "level"),
    [
        (79.99, Level.NOT_COMPLIANT),
        (80.0, Level.WARNING),
        (95.0, Level.WARNING),
        (95.01, Level.OK),
    ],
)
def test_alignment_bounds_copy_the_portal(aligned, level):
    row = {
        "sample_id": "a",
        "alignment_pct": aligned,
        "overlap_peaks": 200000,
        "idr_peaks": 90000,
        "nfr_peak": True,
    }
    assert levels_of([row]) == [level]


@pytest.mark.parametrize(
    ("overlap", "idr", "level"),
    [
        (160000, 10000, Level.OK),  # one passing file suppresses both audits
        (10000, 80000, Level.OK),
        (120000, 10000, Level.WARNING),
        (90000, 60000, Level.WARNING),
        (90000, 40000, Level.NOT_COMPLIANT),
    ],
)
def test_reproducible_peaks_audit_reads_both_files(overlap, idr, level):
    row = {
        "sample_id": "a",
        "alignment_pct": 99.0,
        "overlap_peaks": overlap,
        "idr_peaks": idr,
        "nfr_peak": True,
    }
    assert levels_of([row]) == [level]


def test_null_metric_fires_nothing():
    row = {
        "sample_id": "a",
        "alignment_pct": None,
        "overlap_peaks": 200000,
        "idr_peaks": 90000,
        "nfr_peak": None,
    }
    assert levels_of([row]) == [Level.OK]


def test_cell_takes_the_highest_level_reading_it():
    metrics = pl.DataFrame(
        {
            "sample_id": ["a"],
            "alignment_pct": [70.0],
            "overlap_peaks": [120000],
            "idr_peaks": [10000],
            "nfr_peak": [False],
        }
    )
    cells = cell_levels(run_audits(metrics, atac_spec()))
    assert cells[("a", "alignment_pct")] == Level.NOT_COMPLIANT
    assert cells[("a", "overlap_peaks")] == Level.WARNING
    assert cells[("a", "idr_peaks")] == Level.WARNING
    assert cells[("a", "nfr_peak")] == Level.NOT_COMPLIANT


def test_bounds_fall_on_the_cutoffs():
    spec = atac_spec()
    aligned = spec.metrics[0]
    bounds = audit_bounds(aligned, spec, 60, 100)
    assert bounds == pytest.approx([80, 95], abs=0.05)
    assert audit_bounds(spec.metrics[1], spec, 0, 300000) == []
