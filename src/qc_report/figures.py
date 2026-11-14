"""Plotly figures: cohort histograms and per-sample panels from common tool outputs.

Figures carry no template. The report applies the spec theme's poikilos
template once, in the browser.
"""

import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots

from .spec import Level

# ENCODE's audit colors (green, yellow, orange, red) in zed swatches.
LEVEL_COLORS = {
    Level.OK: "#00a651",
    Level.WARNING: "#f2c80f",
    Level.NOT_COMPLIANT: "#f7941d",
    Level.ERROR: "#ed1c24",
}
MUTED = "#a7a9ac"


def format_count(value):
    """Print a count with three significant digits and a k, M or B suffix.

    Parameters
    ----------
    value : float
        Count.

    Returns
    -------
    str
        The count, e.g. ``"38.3M"``.

    Examples
    --------
    >>> [format_count(v) for v in (950, 285431, 38266372)]
    ['950', '285k', '38.3M']
    """
    for size, suffix in ((1e9, "B"), (1e6, "M"), (1e3, "k")):
        if abs(value) >= size:
            return f"{value / size:.3g}{suffix}"
    return f"{value:.0f}"


def level_histogram(values, levels, bounds=(), log=False):
    """Draw a cohort histogram with each bar stacked by audit level.

    Bins are counted here, so the page gets counts, not every value. A log
    histogram bins log10 values and labels its axis in counts.

    Parameters
    ----------
    values : sequence of float
        One value per sample. Non-finite values are dropped.
    levels : sequence of Level
        Each sample's level for this metric.
    bounds : sequence of float, optional
        Audit cutoffs, drawn as dotted lines.
    log : bool, default False
        Bin on log10 values. Non-positive values are dropped.

    Returns
    -------
    plotly.graph_objects.Figure
        Stacked bars, one trace per level present.

    Examples
    --------
    >>> levels = [Level.OK, Level.OK, Level.WARNING, Level.OK]
    >>> fig = level_histogram([9.0, 8.5, 6.0, 12.0], levels, bounds=[7.0])
    >>> [(trace.name, int(sum(trace.y))) for trace in fig.data]
    [('OK', 3), ('WARNING', 1)]
    """
    x = np.asarray(values, dtype=float)
    level = np.asarray([int(v) for v in levels])
    keep = np.isfinite(x) & (x > 0 if log else True)
    x, level = x[keep], level[keep]
    fig = go.Figure()
    if x.size == 0:
        return fig
    scaled = np.log10(x) if log else x
    edges = np.histogram_bin_edges(scaled, bins=min(40, max(10, x.size // 10)))
    ranges = np.column_stack([edges[:-1], edges[1:]])
    if log:
        ranges = 10**ranges
    number = ".3s" if log else ".3g"
    for value in sorted(set(level)):
        counts, _ = np.histogram(scaled[level == value], bins=edges)
        fig.add_bar(
            x=(edges[:-1] + edges[1:]) / 2,
            y=counts,
            width=np.diff(edges),
            customdata=ranges,
            marker={"color": LEVEL_COLORS[Level(value)], "line": {"width": 0}},
            name=Level(value).label,
            hovertemplate=f"%{{customdata[0]:{number}}} to %{{customdata[1]:{number}}}"
            "<br>%{y} samples<extra>%{fullData.name}</extra>",
        )
    for bound in bounds:
        position = np.log10(bound) if log else bound
        fig.add_vline(x=position, line={"dash": "dot", "width": 1, "color": MUTED})
    fig.update_layout(barmode="stack", bargap=0, showlegend=False)
    fig.update_yaxes(title_text="Samples")
    if log:
        ticks = _log_ticks(edges[0], edges[-1])
        fig.update_xaxes(
            tickvals=np.log10(ticks), ticktext=[format_count(t) for t in ticks]
        )
    return fig


def _log_ticks(low, high):
    # 1-2-5 steps between two log10 values; add 3 and 7 when the range is narrow.
    decades = range(int(np.floor(low)), int(np.ceil(high)) + 1)
    for steps in ((1,), (1, 2, 5), (1, 2, 3, 5, 7)):
        ticks = [
            s * 10.0**d
            for d in decades
            for s in steps
            if low <= np.log10(s * 10.0**d) <= high
        ]
        if len(ticks) >= 3:
            return ticks
    return ticks


def read_quality_figure(fastp_json):
    """Draw mean base quality per cycle for each mate, from a fastp report.

    Parameters
    ----------
    fastp_json : str or pathlib.Path
        Report from ``fastp --json``.

    Returns
    -------
    plotly.graph_objects.Figure
        One line per mate, before filtering.

    Examples
    --------
    >>> import json
    >>> curves = {"quality_curves": {"mean": [32.1, 34.0, 33.5]}}
    >>> report = {"read1_before_filtering": curves, "read2_before_filtering": curves}
    >>> path = tmp_path / "fastp.json"
    >>> _ = path.write_text(json.dumps(report))
    >>> [trace.name for trace in read_quality_figure(path).data]
    ['Read 1', 'Read 2']
    """
    report = json.loads(Path(fastp_json).read_text())
    fig = go.Figure()
    for mate in ("read1", "read2"):
        curves = report.get(f"{mate}_before_filtering", {}).get("quality_curves")
        if curves:
            mean = curves["mean"]
            fig.add_scatter(
                x=np.arange(1, len(mean) + 1),
                y=mean,
                mode="lines",
                name=mate.replace("read", "Read "),
            )
    fig.update_layout(xaxis_title="Cycle", yaxis_title="Mean quality")
    return fig


def gc_bias_figure(gc_metrics):
    """Draw Picard's normalized coverage by GC over the genome's GC windows.

    Parameters
    ----------
    gc_metrics : polars.DataFrame
        CollectGcBiasMetrics detail table, read with `read_picard_metrics`.

    Returns
    -------
    plotly.graph_objects.Figure
        Coverage with error bars, the genome's windows as bars behind it, and
        a dotted line at even coverage.

    Examples
    --------
    >>> import polars as pl
    >>> gc = pl.DataFrame(
    ...     {
    ...         "GC": [40, 50, 60],
    ...         "WINDOWS": [300, 500, 200],
    ...         "NORMALIZED_COVERAGE": [0.9, 1.1, 1.3],
    ...         "ERROR_BAR_WIDTH": [0.01, 0.01, 0.02],
    ...     }
    ... )
    >>> [trace.name for trace in gc_bias_figure(gc).data]
    ['Genome windows', 'Normalized coverage']
    """
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    windows = gc_metrics["WINDOWS"] / gc_metrics["WINDOWS"].sum()
    fig.add_bar(
        x=gc_metrics["GC"],
        y=windows,
        name="Genome windows",
        marker_color="#e6e7e8",
        secondary_y=True,
        hovertemplate="GC %{x}%<br>%{y:.2%} of windows<extra></extra>",
    )
    fig.add_scatter(
        x=gc_metrics["GC"],
        y=gc_metrics["NORMALIZED_COVERAGE"],
        error_y={"array": gc_metrics["ERROR_BAR_WIDTH"], "thickness": 0.5, "width": 0},
        mode="lines+markers",
        name="Normalized coverage",
        hovertemplate="GC %{x}%<br>coverage %{y:.2f}<extra></extra>",
    )
    fig.add_hline(y=1, line={"dash": "dot", "width": 1, "color": MUTED})
    fig.update_xaxes(title_text="GC (%)")
    fig.update_yaxes(title_text="Normalized coverage", secondary_y=False)
    # The windows are context, so their axis stays bare.
    fig.update_yaxes(
        secondary_y=True, showline=False, showticklabels=False, ticks="", title_text=""
    )
    fig.update_layout(showlegend=False)
    return fig


def tss_profile_figure(profile):
    """Draw the mean signal around TSSs, normalized to the flanks.

    Parameters
    ----------
    profile : polars.DataFrame
        ``distance`` from the TSS in bp, and ``enrichment``, the mean signal
        over the mean of the outer flanks.

    Returns
    -------
    plotly.graph_objects.Figure
        The profile and a dotted line at the TSS. The title gives the peak,
        which is the TSS enrichment score.

    Examples
    --------
    >>> import polars as pl
    >>> profile = pl.DataFrame({"distance": [-10, 0, 10], "enrichment": [1.0, 8.6, 1.1]})
    >>> tss_profile_figure(profile).layout.title.text
    'TSS enrichment 8.60'
    """
    fig = go.Figure(
        go.Scatter(
            x=profile["distance"],
            y=profile["enrichment"].round(4),
            mode="lines",
            hovertemplate="%{x:.0f} bp<br>%{y:.2f}<extra></extra>",
        )
    )
    fig.add_vline(x=0, line={"dash": "dot", "width": 1, "color": MUTED})
    fig.update_layout(
        title_text=f"TSS enrichment {profile['enrichment'].max():.2f}",
        xaxis_title="Distance from TSS (bp)",
        yaxis_title="Enrichment",
        showlegend=False,
    )
    return fig


def insert_size_figure(histogram, bands=(), longest=1000):
    """Draw fragment counts per length from Picard's insert size histogram.

    Parameters
    ----------
    histogram : polars.DataFrame
        The ``## HISTOGRAM`` table of CollectInsertSizeMetrics, read with
        ``read_picard_metrics(path, histogram=True)``. The first count column
        is drawn.
    bands : sequence of tuple, optional
        ``(label, low, high)`` ranges in bp to shade, such as the
        nucleosome-free range.
    longest : int, default 1000
        Longest fragment shown.

    Returns
    -------
    plotly.graph_objects.Figure
        Filled step line with the bands behind it.

    Examples
    --------
    >>> import polars as pl
    >>> histogram = pl.DataFrame(
    ...     {"insert_size": [50, 51, 200], "All_Reads.fr_count": [7, 9, 30]}
    ... )
    >>> fig = insert_size_figure(histogram, bands=[("NFR", 0, 150)])
    >>> [shape.x1 for shape in fig.layout.shapes]
    [150]
    """
    sizes = histogram["insert_size"]
    counts = histogram[histogram.columns[1]]
    fig = go.Figure(
        go.Scatter(
            x=sizes,
            y=counts,
            mode="lines",
            line_shape="hvh",
            fill="tozeroy",
            hovertemplate="%{x} bp<br>%{y:,} fragments<extra></extra>",
        )
    )
    for index, (label, low, high) in enumerate(bands):
        # Alternate shades, so touching bands stay apart.
        fig.add_vrect(
            x0=low,
            x1=high,
            fillcolor=MUTED,
            opacity=0.12 if index % 2 == 0 else 0.24,
            line_width=0,
            layer="below",
            annotation_text=label,
            annotation_position="top left",
        )
    fig.update_xaxes(
        range=[0, min(longest, sizes.max())], title_text="Fragment length (bp)"
    )
    fig.update_yaxes(title_text="Fragments")
    fig.update_layout(showlegend=False)
    return fig


def fingerprint_figure(curve):
    """Draw a fingerprint: the share of reads in the least covered share of bins.

    Parameters
    ----------
    curve : polars.DataFrame
        ``fraction_of_bins`` and ``fraction_of_reads``, both increasing from
        0 to 1, with bins sorted by read count as deepTools plotFingerprint
        sorts them.

    Returns
    -------
    plotly.graph_objects.Figure
        The curve and a dotted diagonal, which uniform coverage would follow.
        Enriched libraries bow further below it.

    Examples
    --------
    >>> import polars as pl
    >>> curve = pl.DataFrame(
    ...     {"fraction_of_bins": [0.0, 0.5, 1.0], "fraction_of_reads": [0.0, 0.2, 1.0]}
    ... )
    >>> len(fingerprint_figure(curve).data)
    1
    """
    fig = go.Figure(
        go.Scatter(
            x=curve["fraction_of_bins"].round(4),
            y=curve["fraction_of_reads"].round(4),
            mode="lines",
            hovertemplate="%{x:.0%} of bins<br>%{y:.1%} of reads<extra></extra>",
        )
    )
    fig.add_shape(
        type="line",
        x0=0,
        y0=0,
        x1=1,
        y1=1,
        line={"dash": "dot", "width": 1, "color": MUTED},
    )
    fig.update_xaxes(title_text="Fraction of bins", range=[0, 1])
    fig.update_yaxes(title_text="Fraction of reads", range=[0, 1])
    fig.update_layout(showlegend=False)
    return fig


def idr_signal_figure(idr_txt, threshold=0.05, max_points=2500, seed=0):
    r"""Draw each peak's signal in the two pseudoreplicates, split by IDR.

    Reproducible peaks sit on the diagonal and pass. The legend counts every
    peak. The scatter draws a random sample, so a large file stays small.

    Parameters
    ----------
    idr_txt : str or pathlib.Path
        IDR 2.x output, optionally gzipped. Column 12 holds the global
        -log10 IDR. Columns 15 and 19 hold the signal each replicate was
        ranked on. Peaks missing from a replicate are skipped.
    threshold : float, default 0.05
        IDR threshold.
    max_points : int, default 2500
        Peaks drawn.
    seed : int, default 0
        Sampling seed, so a rerun draws the same peaks.

    Returns
    -------
    plotly.graph_objects.Figure
        Log-log scatter: failing peaks muted, passing peaks in the theme's
        first color, drawn on top.

    Examples
    --------
    >>> def peak(start, idr, signal1, signal2):
    ...     fields = ["chr1", start, start + 500, ".", 1000, ".", 9.0, 50.0, 45.0, 250]
    ...     fields += [idr, idr, start, start + 500, signal1, 250]
    ...     fields += [start, start + 500, signal2, 250]
    ...     return "\t".join(map(str, fields))
    >>> path = tmp_path / "idr.txt"
    >>> _ = path.write_text(
    ...     "\n".join([peak(100, 5.0, 90.0, 80.0), peak(900, 3.2, 40.0, 50.0),
    ...                peak(2000, 0.4, 10.0, 12.0)]) + "\n"
    ... )
    >>> [trace.name for trace in idr_signal_figure(path).data]
    ['Fail IDR 0.05: 1', 'Pass IDR 0.05: 2']
    """
    table = pl.read_csv(idr_txt, separator="\t", has_header=False)
    peaks = table.select(
        passes=pl.nth(11) >= -np.log10(threshold),
        rep1=pl.nth(14),
        rep2=pl.nth(18),
    ).filter((pl.col("rep1") > 0) & (pl.col("rep2") > 0))
    shown = peaks.sample(max_points, seed=seed) if peaks.height > max_points else peaks
    fig = go.Figure()
    for passes, word, marker in ((False, "Fail", {"color": MUTED}), (True, "Pass", {})):
        part = shown.filter(pl.col("passes") == passes)
        total = peaks.filter(pl.col("passes") == passes).height
        fig.add_scatter(
            x=part["rep1"].round(2),
            y=part["rep2"].round(2),
            mode="markers",
            name=f"{word} IDR {threshold:g}: {total:,}",
            marker={"size": 3, "line": {"width": 0}, **marker},
            hovertemplate="%{x:.1f}, %{y:.1f}<extra></extra>",
        )
    fig.update_xaxes(type="log", dtick=1, title_text="Signal, pseudoreplicate 1")
    fig.update_yaxes(type="log", dtick=1, title_text="Signal, pseudoreplicate 2")
    fig.update_layout(
        legend={"x": 0.02, "y": 0.98, "xanchor": "left", "yanchor": "top"}
    )
    return fig
