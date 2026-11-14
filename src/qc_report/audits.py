"""Run a spec's audits over a metrics table."""

import polars as pl

from .spec import Level

FLAG_SCHEMA = {
    "sample_id": pl.String,
    "flag": pl.String,
    "level": pl.Int8,
    "when": pl.String,
    "metrics": pl.List(pl.String),
}


def run_audits(metrics, spec):
    """Find the audits that fire on each sample.

    Parameters
    ----------
    metrics : polars.DataFrame
        One row per sample: ``sample_id`` and a column per metric id.
    spec : Spec
        Audits to run.

    Returns
    -------
    polars.DataFrame
        One row per fired audit: ``sample_id``, ``flag``, ``level`` (a
        `Level` value), ``when`` and the ``metrics`` the condition reads.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Audit, Group, Metric, Spec
    >>> spec = Spec(
    ...     "QC",
    ...     [Group("Enrichment", [Metric("frip", "FRiP")])],
    ...     [Audit("low FRiP score", "NOT_COMPLIANT", "frip < 0.2")],
    ... )
    >>> metrics = pl.DataFrame({"sample_id": ["a", "b"], "frip": [0.1, 0.4]})
    >>> run_audits(metrics, spec)["sample_id"].to_list()
    ['a']
    """
    fired = [
        metrics.filter(audit.condition().fill_null(False))
        .select("sample_id")
        .with_columns(
            flag=pl.lit(audit.flag),
            level=pl.lit(int(audit.severity), dtype=pl.Int8),
            when=pl.lit(audit.when),
            metrics=pl.lit(audit.columns(), dtype=pl.List(pl.String)),
        )
        for audit in spec.audits
    ]
    return pl.concat(fired) if fired else pl.DataFrame(schema=FLAG_SCHEMA)


def cell_levels(flags):
    """Give each metric of each sample the highest level among its fired audits.

    Parameters
    ----------
    flags : polars.DataFrame
        Output of `run_audits`.

    Returns
    -------
    dict of (str, str) to Level
        Level per ``(sample_id, metric)``. Pairs without a fired audit are
        absent, which means `Level.OK`.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Audit, Group, Metric, Spec
    >>> spec = Spec(
    ...     "QC",
    ...     [Group("Library", [Metric("pbc1", "PBC1")])],
    ...     [
    ...         Audit("severe bottlenecking", "NOT_COMPLIANT", "pbc1 < 0.7"),
    ...         Audit("mild to moderate bottlenecking", "WARNING", "pbc1 <= 0.9"),
    ...     ],
    ... )
    >>> flags = run_audits(pl.DataFrame({"sample_id": ["a"], "pbc1": [0.5]}), spec)
    >>> cell_levels(flags)[("a", "pbc1")].name
    'NOT_COMPLIANT'
    """
    levels = (
        flags.explode("metrics")
        .group_by("sample_id", "metrics")
        .agg(pl.col("level").max())
        .iter_rows()
    )
    return {(sample, metric): Level(level) for sample, metric, level in levels}


def sample_levels(flags):
    """Give each sample the highest level among its fired audits, as the portal does.

    Parameters
    ----------
    flags : polars.DataFrame
        Output of `run_audits`.

    Returns
    -------
    dict of str to Level
        Level per sample with a fired audit. Absent samples are `Level.OK`.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Audit, Group, Metric, Spec
    >>> spec = Spec(
    ...     "QC",
    ...     [Group("Enrichment", [Metric("tss", "TSS"), Metric("frip", "FRiP")])],
    ...     [
    ...         Audit("low TSS enrichment", "NOT_COMPLIANT", "tss < 5"),
    ...         Audit("moderate FRiP score", "WARNING", "frip <= 0.3"),
    ...     ],
    ... )
    >>> metrics = pl.DataFrame(
    ...     {"sample_id": ["a", "b", "c"], "tss": [4.0, 9.0, 9.0], "frip": [0.2, 0.2, 0.5]}
    ... )
    >>> {s: level.name for s, level in sorted(sample_levels(run_audits(metrics, spec)).items())}
    {'a': 'NOT_COMPLIANT', 'b': 'WARNING'}
    """
    levels = flags.group_by("sample_id").agg(pl.col("level").max()).iter_rows()
    return {sample: Level(level) for sample, level in levels}


def audit_bounds(metric, spec, low, high, points=2001):
    """Find where the audits of one metric change level between two values.

    Only audits that read this metric alone count, so the bounds can be
    drawn on its histogram.

    Parameters
    ----------
    metric : Metric
        Metric to scan.
    spec : Spec
        Audits to evaluate.
    low, high : float
        Range to scan.
    points : int, default 2001
        Grid size. A bound is placed within ``(high - low) / (points - 1)``.

    Returns
    -------
    list of float
        Values where the level changes, in increasing order.

    Examples
    --------
    >>> from qc_report.spec import Audit, Group, Metric, Spec
    >>> tss = Metric("tss", "TSS")
    >>> spec = Spec(
    ...     "QC",
    ...     [Group("Enrichment", [tss])],
    ...     [
    ...         Audit("low TSS enrichment", "NOT_COMPLIANT", "tss < 5"),
    ...         Audit("moderate TSS enrichment", "WARNING", "tss >= 5 AND tss <= 7"),
    ...     ],
    ... )
    >>> [round(bound, 1) for bound in audit_bounds(tss, spec, 0, 20)]
    [5.0, 7.0]
    """
    own = [audit for audit in spec.audits if audit.columns() == [metric.id]]
    if not own or not low < high:
        return []
    grid = pl.DataFrame({metric.id: pl.linear_space(low, high, points, eager=True)})
    level = pl.max_horizontal(
        pl.when(audit.condition()).then(int(audit.severity)).otherwise(0)
        for audit in own
    )
    scanned = grid.with_columns(level=level)
    changes = scanned.filter(pl.col("level").diff() != 0)
    step = (high - low) / (points - 1)
    return [value - step / 2 for value in changes[metric.id]]
