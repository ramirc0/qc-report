"""Read common tool outputs into flat metric dicts and tables."""

import io
import json
from pathlib import Path

import polars as pl
from zlib_ng import gzip_ng


def read_fastp(path):
    """Read raw-read metrics from a fastp JSON report.

    Parameters
    ----------
    path : str or pathlib.Path
        Report from ``fastp --json``.

    Returns
    -------
    dict
        ``reads`` (both mates counted), ``read_length`` (mate 1 mean, before
        filtering), ``q30_pct``, ``gc_pct`` and ``duplication_pct``.

    Examples
    --------
    >>> import json
    >>> before = {
    ...     "total_reads": 2000,
    ...     "read1_mean_length": 76,
    ...     "q30_rate": 0.93,
    ...     "gc_content": 0.47,
    ... }
    >>> report = {"summary": {"before_filtering": before}, "duplication": {"rate": 0.02}}
    >>> path = tmp_path / "fastp.json"
    >>> _ = path.write_text(json.dumps(report))
    >>> read_fastp(path)["reads"]
    2000
    """
    report = json.loads(Path(path).read_text())
    before = report["summary"]["before_filtering"]
    return {
        "reads": before["total_reads"],
        "read_length": before["read1_mean_length"],
        "q30_pct": 100 * before["q30_rate"],
        "gc_pct": 100 * before["gc_content"],
        "duplication_pct": 100 * report["duplication"]["rate"],
    }


def read_cutadapt(path):
    """Read trimming metrics from a cutadapt JSON report.

    Parameters
    ----------
    path : str or pathlib.Path
        Report from ``cutadapt --json``.

    Returns
    -------
    dict
        ``reads_in`` and ``reads_out`` (pairs for paired-end input),
        ``adapter_pct`` (mate 1 reads with an adapter) and ``bases_kept_pct``.

    Examples
    --------
    >>> import json
    >>> report = {
    ...     "read_counts": {"input": 100, "output": 98, "read1_with_adapter": 25},
    ...     "basepair_counts": {"input": 7600, "output": 7000},
    ... }
    >>> path = tmp_path / "cutadapt.json"
    >>> _ = path.write_text(json.dumps(report))
    >>> read_cutadapt(path)["adapter_pct"]
    25.0
    """
    report = json.loads(Path(path).read_text())
    reads = report["read_counts"]
    bases = report["basepair_counts"]
    return {
        "reads_in": reads["input"],
        "reads_out": reads["output"],
        "adapter_pct": 100 * reads["read1_with_adapter"] / reads["input"],
        "bases_kept_pct": 100 * bases["output"] / bases["input"],
    }


def read_flagstat(path):
    """Read QC-passed counts from ``samtools flagstat --output-fmt json``.

    Parameters
    ----------
    path : str or pathlib.Path
        The flagstat JSON.

    Returns
    -------
    dict
        ``primary`` reads, ``primary_mapped_pct`` (each read once, mapped if
        any alignment is), ``properly_paired_pct``, ``secondary`` and
        ``duplicates``.

    Examples
    --------
    >>> import json
    >>> passed = {
    ...     "primary": 1000,
    ...     "primary mapped %": 97.5,
    ...     "properly paired %": 95.1,
    ...     "secondary": 40,
    ...     "primary duplicates": 0,
    ... }
    >>> path = tmp_path / "flagstat.json"
    >>> _ = path.write_text(json.dumps({"QC-passed reads": passed}))
    >>> read_flagstat(path)["primary_mapped_pct"]
    97.5
    """
    passed = json.loads(Path(path).read_text())["QC-passed reads"]
    return {
        "primary": passed["primary"],
        "primary_mapped_pct": passed["primary mapped %"],
        "properly_paired_pct": passed["properly paired %"],
        "secondary": passed["secondary"],
        "duplicates": passed["primary duplicates"],
    }


def read_picard_metrics(path, histogram=False):
    r"""Read a Picard metrics file's metrics table or its histogram.

    Parameters
    ----------
    path : str or pathlib.Path
        Picard output, such as CollectGcBiasMetrics or
        CollectInsertSizeMetrics.
    histogram : bool, default False
        Read the ``## HISTOGRAM`` table instead of the ``## METRICS CLASS``
        table.

    Returns
    -------
    polars.DataFrame
        The table, with Picard's column names.

    Raises
    ------
    ValueError
        If the file has no such table.

    Examples
    --------
    >>> path = tmp_path / "insert_size.txt"
    >>> _ = path.write_text(
    ...     "## htsjdk.samtools.metrics.StringHeader\n\n"
    ...     "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n"
    ...     "MEDIAN_INSERT_SIZE\tMODE_INSERT_SIZE\n120\t60\n\n"
    ...     "## HISTOGRAM\tjava.lang.Integer\n"
    ...     "insert_size\tAll_Reads.fr_count\n50\t7\n51\t9\n"
    ... )
    >>> read_picard_metrics(path)["MEDIAN_INSERT_SIZE"].item()
    120
    >>> read_picard_metrics(path, histogram=True).height
    2
    """
    marker = "## HISTOGRAM" if histogram else "## METRICS CLASS"
    lines = Path(path).read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith(marker)]
    if not starts:
        raise ValueError(f"{path}: no {marker!r} table")
    table = []
    for line in lines[starts[0] + 1 :]:
        if not line.strip() or line.startswith("#"):
            break
        table.append(line)
    return pl.read_csv(io.StringIO("\n".join(table)), separator="\t")


def read_fingerprint(path):
    r"""Read deepTools plotFingerprint ``--outQualityMetrics`` for one sample.

    Parameters
    ----------
    path : str or pathlib.Path
        Metrics table with one row.

    Returns
    -------
    dict
        ``auc``, ``synthetic_auc``, ``x_intercept``, ``elbow_point`` and
        ``synthetic_jsd``.

    Examples
    --------
    >>> path = tmp_path / "fingerprint.tsv"
    >>> _ = path.write_text(
    ...     "Sample\tAUC\tSynthetic AUC\tX-intercept\tSynthetic X-intercept\t"
    ...     "Elbow Point\tSynthetic Elbow Point\tSynthetic JS Distance\n"
    ...     "s1\t0.23\t0.49\t0.18\t0.0\t0.61\t0.51\t0.31\n"
    ... )
    >>> read_fingerprint(path)["synthetic_jsd"]
    0.31
    """
    row = pl.read_csv(path, separator="\t").row(0, named=True)
    return {
        "auc": row["AUC"],
        "synthetic_auc": row["Synthetic AUC"],
        "x_intercept": row["X-intercept"],
        "elbow_point": row["Elbow Point"],
        "synthetic_jsd": row["Synthetic JS Distance"],
    }


def count_records(path):
    r"""Count the lines of a text file, gzipped or not.

    Parameters
    ----------
    path : str or pathlib.Path
        File such as a ``.narrowPeak.gz``.

    Returns
    -------
    int
        Lines in the file.

    Examples
    --------
    >>> from zlib_ng import gzip_ng
    >>> path = tmp_path / "peaks.narrowPeak.gz"
    >>> with gzip_ng.open(path, "wt") as f:
    ...     _ = f.write("chr1\t10\t60\nchr1\t90\t150\n")
    >>> count_records(path)
    2
    """
    path = Path(path)
    opener = gzip_ng.open if path.suffix == ".gz" else open
    with opener(path, "rb") as f:
        return sum(chunk.count(b"\n") for chunk in iter(lambda: f.read(1 << 20), b""))
