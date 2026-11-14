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


def level_histogram(values, levels, bounds=(), log=False):
    """Draw a cohort histogram with each bar stacked by audit level.

    Parameters
    ----------
    values : sequence of float
        One value per sample. Non-finite values are dropped.
    levels : sequence of Level
        Each sample's level for this metric.
    bounds : sequence of float, optional
        Audit cutoffs, drawn as dotted lines.
    log : bool, default False
        Bin and draw on a log10 axis. Non-positive values are dropped.

    Returns
    -------
    plotly.graph_objects.Figure
        Stacked histogram, one trace per level present.

    Examples
    --------
    >>> levels = [Level.OK, Level.OK, Level.WARNING, Level.OK]
    >>> fig = level_histogram([9.0, 8.5, 6.0, 12.0], levels, bounds=[7.0])
    >>> [trace.name for trace in fig.data]
    ['OK', 'WARNING']
    """
    x = np.asarray(values, dtype=float)
    level = np.asarray([int(v) for v in levels])
    keep = np.isfinite(x) & (x > 0 if log else True)
    x, level = x[keep], level[keep]
    fig = go.Figure()
    if x.size == 0:
        return fig
    binned = np.log10(x) if log else x
    edges = np.histogram_bin_edges(binned, bins=min(40, max(10, x.size // 10)))
    if log:
        edges = 10**edges
    for value in sorted(set(level)):
        fig.add_histogram(
            x=x[level == value],
            xbins={"start": edges[0], "end": edges[-1], "size": _bin_size(edges, log)},
            marker_color=LEVEL_COLORS[Level(value)],
            name=Level(value).label,
            hovertemplate="%{x}<br>%{y} samples<extra>%{fullData.name}</extra>",
        )
    for bound in bounds:
        fig.add_vline(x=bound, line={"dash": "dot", "width": 1, "color": MUTED})
    fig.update_layout(barmode="stack", bargap=0.08, showlegend=False)
    fig.update_xaxes(type="log" if log else "linear")
    fig.update_yaxes(title_text="Samples")
    return fig


def _bin_size(edges, log):
    # Plotly sizes log-axis bins in log10 units.
    return float(
        np.log10(edges[1]) - np.log10(edges[0]) if log else edges[1] - edges[0]
    )


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
    fig.update_yaxes(
        title_text="Genome windows", tickformat=".0%", secondary_y=True, showline=False
    )
    fig.update_layout(showlegend=False)
    return fig


def idr_rank_figure(idr_txt, threshold=0.05, bins=100):
    r"""Draw how two pseudoreplicates rank the same peaks, as a 2D histogram.

    Parameters
    ----------
    idr_txt : str or pathlib.Path
        IDR 2.x output, optionally gzipped. Column 12 holds the global
        -log10 IDR; columns 15 and 19 hold each replicate's signal.
    threshold : float, default 0.05
        IDR threshold, for the count in the title.
    bins : int, default 100
        Bins per axis.

    Returns
    -------
    plotly.graph_objects.Figure
        Heatmap of peak counts by rank in each replicate, rank 1 the
        strongest. Color is log10(count + 1).

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
    >>> idr_rank_figure(path, threshold=0.05, bins=3).layout.title.text
    '2 of 3 peaks pass IDR 0.05'
    """
    table = pl.read_csv(idr_txt, separator="\t", has_header=False)
    peaks = table.select(
        idr=pl.nth(11),
        rank1=pl.nth(14).rank("ordinal", descending=True),
        rank2=pl.nth(18).rank("ordinal", descending=True),
    )
    counts, xedges, yedges = np.histogram2d(peaks["rank1"], peaks["rank2"], bins=bins)
    passed = peaks.filter(pl.col("idr") >= -np.log10(threshold)).height
    fig = go.Figure(
        go.Heatmap(
            x=(xedges[:-1] + xedges[1:]) / 2,
            y=(yedges[:-1] + yedges[1:]) / 2,
            z=np.log10(counts.T + 1),
            customdata=counts.T,
            colorscale="Greys",
            showscale=False,
            hovertemplate="rank %{x:.0f} vs %{y:.0f}<br>%{customdata:.0f} peaks"
            "<extra></extra>",
        )
    )
    fig.update_layout(
        title_text=f"{passed:,} of {peaks.height:,} peaks pass IDR {threshold:g}",
        xaxis_title="Rank, pseudoreplicate 1",
        yaxis_title="Rank, pseudoreplicate 2",
    )
    return fig
