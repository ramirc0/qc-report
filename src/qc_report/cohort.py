"""Place each sample within its run: cohort outliers and percentiles."""

import polars as pl

# Iglewicz and Hoaglin's cutoff for the modified z-score.
OUTLIER_Z = 3.5
# Below this many samples, a median absolute deviation says little.
MIN_SAMPLES = 10

OUTLIER_SCHEMA = {
    "sample_id": pl.String,
    "metric": pl.String,
    "side": pl.String,
    "z": pl.Float64,
}


def find_outliers(metrics, spec):
    """Flag samples far from the run's median, on the side a metric marks bad.

    The modified z-score is ``0.6745 * (x - median) / MAD``. Metrics with
    ``scale: log`` are scored on log10 values, and non-positive values are
    skipped. A metric with a MAD of 0 flags nothing.

    Parameters
    ----------
    metrics : polars.DataFrame
        One row per sample: ``sample_id`` and a column per metric id.
    spec : Spec
        Metrics, with the side each flags.

    Returns
    -------
    polars.DataFrame
        One row per outlier: ``sample_id``, ``metric``, ``side`` (``"low"`` or
        ``"high"``) and ``z``. Empty below `MIN_SAMPLES` samples.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Group, Metric, Spec
    >>> spec = Spec("QC", [Group("Depth", [Metric("reads", "Reads", outlier="low")])])
    >>> reads = [100.0 + i for i in range(19)] + [10.0]
    >>> metrics = pl.DataFrame({"sample_id": [f"s{i}" for i in range(20)], "reads": reads})
    >>> find_outliers(metrics, spec)["sample_id"].to_list()
    ['s19']
    """
    if metrics.height < MIN_SAMPLES:
        return pl.DataFrame(schema=OUTLIER_SCHEMA)
    frames = []
    for metric in spec.metrics:
        if metric.outlier == "none" or metric.format == "flag":
            continue
        x = pl.col(metric.id).cast(pl.Float64)
        if metric.scale == "log":
            x = pl.when(x > 0).then(x.log10())
        deviation = x - x.median()
        z = 0.6745 * deviation / deviation.abs().median()
        scored = metrics.select("sample_id", z=z).filter(pl.col("z").is_finite())
        sides = {
            "low": pl.col("z") < -OUTLIER_Z,
            "high": pl.col("z") > OUTLIER_Z,
        }
        for side, test in sides.items():
            if metric.outlier in (side, "both"):
                frames.append(
                    scored.filter(test).select(
                        "sample_id",
                        metric=pl.lit(metric.id),
                        side=pl.lit(side),
                        z="z",
                    )
                )
    return pl.concat(frames) if frames else pl.DataFrame(schema=OUTLIER_SCHEMA)


def percentiles(metrics, spec):
    """Rank each sample within the run, per numeric metric.

    Parameters
    ----------
    metrics : polars.DataFrame
        One row per sample: ``sample_id`` and a column per metric id.
    spec : Spec
        Metrics to rank. Flags are skipped.

    Returns
    -------
    polars.DataFrame
        ``sample_id`` and one column per numeric metric, holding the
        percentile from 0 to 100. Ties share their average rank.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Group, Metric, Spec
    >>> spec = Spec("QC", [Group("Enrichment", [Metric("tss", "TSS")])])
    >>> metrics = pl.DataFrame({"sample_id": ["a", "b"], "tss": [12.0, 6.0]})
    >>> percentiles(metrics, spec)["tss"].to_list()
    [100.0, 50.0]
    """
    numeric = [metric.id for metric in spec.metrics if metric.format != "flag"]
    return metrics.select(
        "sample_id",
        *(
            (100 * pl.col(m).rank("average") / pl.col(m).count()).alias(m)
            for m in numeric
        ),
    )
