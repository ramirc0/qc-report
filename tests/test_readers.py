import json

import pytest
from zlib_ng import gzip_ng

from qc_report import (
    count_records,
    read_cutadapt,
    read_fastp,
    read_fingerprint,
    read_flagstat,
    read_picard_metrics,
)


def test_fastp(tmp_path):
    report = {
        "summary": {
            "before_filtering": {
                "total_reads": 200,
                "read1_mean_length": 76,
                "q30_rate": 0.93,
                "gc_content": 0.47,
            }
        },
        "duplication": {"rate": 0.02},
    }
    path = tmp_path / "fastp.json"
    path.write_text(json.dumps(report))
    assert read_fastp(path) == {
        "reads": 200,
        "read_length": 76,
        "q30_pct": pytest.approx(93),
        "gc_pct": pytest.approx(47),
        "duplication_pct": pytest.approx(2),
    }


def test_cutadapt(tmp_path):
    report = {
        "read_counts": {"input": 100, "output": 98, "read1_with_adapter": 25},
        "basepair_counts": {"input": 1000, "output": 900},
    }
    path = tmp_path / "cutadapt.json"
    path.write_text(json.dumps(report))
    assert read_cutadapt(path) == {
        "reads_in": 100,
        "reads_out": 98,
        "adapter_pct": 25,
        "bases_kept_pct": 90,
    }


def test_flagstat(tmp_path):
    passed = {
        "primary": 1000,
        "primary mapped %": 97.5,
        "properly paired %": 95.1,
        "secondary": 40,
        "primary duplicates": 0,
    }
    path = tmp_path / "flagstat.json"
    path.write_text(json.dumps({"QC-passed reads": passed, "QC-failed reads": {}}))
    assert read_flagstat(path)["primary_mapped_pct"] == 97.5


def test_picard_metrics_and_histogram(tmp_path):
    path = tmp_path / "insert.txt"
    path.write_text(
        "## htsjdk.samtools.metrics.StringHeader\n# CollectInsertSizeMetrics\n\n"
        "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n"
        "MEDIAN_INSERT_SIZE\tMODE_INSERT_SIZE\n120\t60\n\n"
        "## HISTOGRAM\tjava.lang.Integer\ninsert_size\tAll_Reads.fr_count\n"
        "50\t7\n51\t9\n"
    )
    assert read_picard_metrics(path).row(0) == (120, 60)
    assert read_picard_metrics(path, histogram=True)["All_Reads.fr_count"].sum() == 16


def test_fingerprint(tmp_path):
    path = tmp_path / "jsd.qc"
    path.write_text(
        "Sample\tAUC\tSynthetic AUC\tX-intercept\tSynthetic X-intercept\t"
        "Elbow Point\tSynthetic Elbow Point\tSynthetic JS Distance\n"
        "s\t0.23\t0.49\t0.18\t0.0\t0.61\t0.50\t0.31\n"
    )
    assert read_fingerprint(path)["synthetic_jsd"] == 0.31


def test_count_records(tmp_path):
    plain = tmp_path / "peaks.bed"
    plain.write_text("a\nb\nc\n")
    packed = tmp_path / "peaks.bed.gz"
    with gzip_ng.open(packed, "wt") as f:
        f.write("a\nb\n")
    assert count_records(plain) == 3
    assert count_records(packed) == 2
