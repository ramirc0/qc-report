"""Assemble the report page and the per-sample figure scripts."""

import datetime as dt
import math
import re
from importlib import resources
from pathlib import Path

import msgspec
import plotly.graph_objects as go
import plotly.io as pio
import poikilos as pk
import polars as pl
from minijinja import Environment
from plotly.offline import get_plotlyjs

from .audits import audit_bounds, cell_levels, run_audits, sample_levels
from .cohort import find_outliers, percentiles
from .figures import level_histogram
from .spec import Level

SAMPLES_DIR = "samples"
_TEMPLATES = resources.files("qc_report") / "templates"


def format_value(value, metric):
    """Print a metric value the way its spec asks.

    Parameters
    ----------
    value : float, int, bool or None
        The value. None and NaN print as ``n/a``.
    metric : Metric
        Format and digits.

    Returns
    -------
    str
        Counts get three significant digits and a k, M or B suffix.

    Examples
    --------
    >>> from qc_report.spec import Metric
    >>> format_value(38266372, Metric("fragments", "Fragments", format="count"))
    '38.3M'
    >>> format_value(97.651, Metric("aligned", "Aligned", format="percent", digits=1))
    '97.7%'
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    if metric.format == "flag":
        return "yes" if value else "no"
    if isinstance(value, float) and math.isinf(value):
        return "∞" if value > 0 else "-∞"
    if metric.format == "count":
        for size, suffix in ((1e9, "B"), (1e6, "M"), (1e3, "k")):
            if abs(value) >= size:
                return f"{value / size:.3g}{suffix}"
        return f"{value:.0f}"
    unit = "%" if metric.format == "percent" else ""
    return f"{value:.{metric.digits}f}{unit}"


def figure_json(fig):
    """Serialize a figure for the report without its template.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Figure to serialize.

    Returns
    -------
    str
        JSON with ``data`` and ``layout``. The page applies the theme.

    Examples
    --------
    >>> import json
    >>> import plotly.graph_objects as go
    >>> fig = go.Figure(go.Bar(y=[3, 1]), layout={"template": "plotly_white"})
    >>> sorted(json.loads(figure_json(fig))["layout"])
    []
    """
    spec = fig.to_plotly_json()
    spec["layout"].pop("template", None)
    return pio.json.to_json_plotly(spec)


def write_sample_panels(path, sample, panels):
    """Write a sample's figures as the script its detail view loads.

    Parameters
    ----------
    path : str or pathlib.Path
        Output, ``<report dir>/samples/<sample>.js``.
    sample : str
        Sample id, as in the metrics table.
    panels : list of tuple
        ``(title, content)`` pairs in display order. Content is a plotly
        figure, or a path to an image relative to the report's directory.

    Examples
    --------
    >>> import plotly.graph_objects as go
    >>> path = tmp_path / "samples" / "s1.js"
    >>> panels = [("Curve", go.Figure(go.Scatter(y=[1, 2]))), ("Plot", "../s1/tss.png")]
    >>> write_sample_panels(path, "s1", panels)
    >>> path.read_text()[:39]
    'qcReport.addSample({"sample": "s1", "pa'
    """
    items = [
        f'{{"title": {_json(title)}, "figure": {figure_json(content)}}}'
        if isinstance(content, go.Figure)
        else f'{{"title": {_json(title)}, "image": {_json(str(content))}}}'
        for title, content in panels
    ]
    payload = f'{{"sample": {_json(sample)}, "panels": [{", ".join(items)}]}}'
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(f"qcReport.addSample({payload});\n")


def render_report(metrics, spec, run_id, version="", run_href=".."):
    """Render the run's report page.

    Parameters
    ----------
    metrics : polars.DataFrame
        One row per sample: ``sample_id`` and a column per metric id.
    spec : Spec
        Layout, audits, links and theme.
    run_id : str
        Run name for the header.
    version : str, optional
        Pipeline version for the header.
    run_href : str, default ".."
        Path from the report's directory to the run directory, for links.

    Returns
    -------
    str
        A self-contained HTML page. Sample details load
        ``samples/<sample>.js`` beside it when they exist.

    Examples
    --------
    >>> import polars as pl
    >>> from qc_report.spec import Audit, Group, Metric, Spec
    >>> spec = Spec(
    ...     "QC",
    ...     [Group("Enrichment", [Metric("frip", "FRiP", table=True)])],
    ...     [Audit("low FRiP score", "NOT_COMPLIANT", "frip < 0.2")],
    ... )
    >>> metrics = pl.DataFrame({"sample_id": ["a", "b"], "frip": [0.1, 0.4]})
    >>> page = render_report(metrics, spec, "run1")
    >>> 'data-sample="a" data-level="2"' in page
    True
    """
    flags = run_audits(metrics, spec)
    cells = cell_levels(flags)
    status = sample_levels(flags)
    outliers = {
        (sample, metric): side
        for sample, metric, side, _ in find_outliers(metrics, spec).iter_rows()
    }
    ranks = {
        row["sample_id"]: row
        for row in percentiles(metrics, spec).iter_rows(named=True)
    }
    fired = _fired_by_sample(flags)
    template = pk.plotly_template(spec.theme)
    rows = list(metrics.iter_rows(named=True))
    sample_view = _SampleView(spec, cells, outliers, ranks, fired, status, run_href)
    cards = [_card(m, metrics, cells, spec) for m in spec.metrics if m.headline]

    figures = ", ".join(f"{_json(card['id'])}: {card.pop('figure')}" for card in cards)
    data_json = (
        f'{{"template": {_json(template.to_plotly_json())}, '
        f'"cards": {{{figures}}}, "samplesDir": {_json(SAMPLES_DIR)}}}'
    )

    context = {
        "title": spec.title,
        "run_id": run_id,
        "version": version,
        "generated": dt.date.today().isoformat(),
        "standards": msgspec.structs.asdict(spec.standards) if spec.standards else None,
        "n_samples": len(rows),
        "font": template.layout.font.family,
        "level_labels": {int(level): level.label for level in Level},
        "status": _status_counts(
            [status.get(r["sample_id"], Level.OK) for r in rows], outliers
        ),
        "flag_summary": _flag_summary(flags, len(rows)),
        "cards": cards,
        "groups": [
            {
                "name": group.name,
                "metrics": [_header(m) for m in group.metrics if m.table],
            }
            for group in spec.groups
            if any(m.table for m in group.metrics)
        ],
        "rows": [sample_view.row(row) for row in rows],
        "details": [sample_view.detail(row) for row in rows],
        "audits": [
            {
                "flag": a.flag,
                "level": int(a.severity),
                "when": a.when,
                "source": a.source,
            }
            for a in spec.audits
        ],
        "metrics": [_header(m) for m in spec.metrics if m.description],
        "data_json": data_json.replace("<", "\\u003c"),
        "css": (_TEMPLATES / "report.css").read_text(),
        "js": (_TEMPLATES / "report.js").read_text(),
        "plotly_js": get_plotlyjs(),
    }
    env = Environment(
        templates={"report.html": (_TEMPLATES / "report.html").read_text()},
        auto_escape_callback=lambda name: True,
    )
    return env.render_template("report.html", **context)


class _SampleView:
    """Table rows and detail views, sharing one run's audit results."""

    def __init__(self, spec, cells, outliers, ranks, fired, status, run_href):
        self.spec = spec
        self.cells = cells
        self.outliers = outliers
        self.ranks = ranks
        self.fired = fired
        self.status = status
        self.run_href = run_href

    def _value(self, row, metric):
        sample = row["sample_id"]
        rank = self.ranks.get(sample, {}).get(metric.id)
        return {
            "text": format_value(row[metric.id], metric),
            "level": int(self.cells.get((sample, metric.id), Level.OK)),
            "outlier": self.outliers.get((sample, metric.id)),
            "percentile": None if rank is None else round(rank),
        }

    def row(self, row):
        sample = row["sample_id"]
        fired = self.fired.get(sample, [])
        per_level = {}
        for flag in fired:
            per_level[flag["level"]] = per_level.get(flag["level"], 0) + 1
        cells = []
        for metric in self.spec.metrics:
            if not metric.table:
                continue
            value = row[metric.id]
            cell = self._value(row, metric)
            reading = [f["flag"] for f in fired if metric.id in f["metrics"]]
            cell["sort"] = "" if value is None else float(value)
            cell["title"] = "; ".join([f"{metric.label} {cell['text']}", *reading])
            cells.append(cell)
        return {
            "sample": sample,
            "level": int(self.status.get(sample, Level.OK)),
            "flags": " ".join(sorted({f["slug"] for f in fired})),
            "flag_names": "; ".join(f["flag"] for f in fired),
            "flag_total": len(fired),
            "flag_counts": [
                {"level": level, "count": n}
                for level, n in sorted(per_level.items(), reverse=True)
            ],
            "outliers": sum(1 for (s, _) in self.outliers if s == sample),
            "cells": cells,
        }

    def detail(self, row):
        sample = row["sample_id"]
        return {
            "sample": sample,
            "level": int(self.status.get(sample, Level.OK)),
            "flags": sorted(self.fired.get(sample, []), key=lambda f: -f["level"]),
            "links": [
                {
                    "label": link.label,
                    "href": f"{self.run_href}/{link.path.format(sample=sample)}",
                }
                for link in self.spec.links
            ],
            "groups": [
                {
                    "name": group.name,
                    "metrics": [
                        {**self._value(row, m), **_header(m)} for m in group.metrics
                    ],
                }
                for group in self.spec.groups
            ],
        }


def _json(value):
    return pio.json.to_json_plotly(value)


def _header(metric):
    return {"id": metric.id, "label": metric.label, "description": metric.description}


def _slug(flag, level):
    return f"{re.sub(r'[^a-z0-9]+', '-', flag.lower()).strip('-')}-{level}"


def _fired_by_sample(flags):
    fired = {}
    for sample, flag, level, when, metrics in flags.iter_rows():
        fired.setdefault(sample, []).append(
            {
                "flag": flag,
                "level": int(level),
                "when": when,
                "metrics": metrics,
                "slug": _slug(flag, level),
            }
        )
    return fired


def _status_counts(levels, outliers):
    total = max(len(levels), 1)
    return {
        "levels": [
            {
                "level": int(level),
                "label": level.label,
                "count": levels.count(level),
                "share": 100 * levels.count(level) / total,
            }
            for level in Level
        ],
        "outliers": len({sample for sample, _ in outliers}),
    }


def _flag_summary(flags, n_samples):
    summary = (
        flags.group_by("flag", "level")
        .agg(
            samples=pl.col("sample_id").n_unique(),
            conditions=pl.col("when").unique(maintain_order=True),
        )
        .sort(["level", "samples"], descending=True)
    )
    return [
        {
            "flag": flag,
            "level": int(level),
            "samples": samples,
            "share": 100 * samples / n_samples,
            "conditions": conditions,
            "slug": _slug(flag, level),
        }
        for flag, level, samples, conditions in summary.iter_rows()
    ]


def _card(metric, metrics, cells, spec):
    values = metrics[metric.id].cast(pl.Float64)
    finite = values.filter(values.is_finite())
    levels = [cells.get((s, metric.id), Level.OK) for s in metrics["sample_id"]]
    bounds = (
        audit_bounds(metric, spec, finite.min(), finite.max())
        if not finite.is_empty()
        else []
    )
    fig = level_histogram(values.to_list(), levels, bounds, log=metric.scale == "log")
    fig.update_xaxes(title_text=metric.label)

    def quantile(q):
        if finite.is_empty():
            return "n/a"
        return format_value(finite.quantile(q), metric)

    return {
        "id": metric.id,
        "label": metric.label,
        "median": quantile(0.5),
        "q1": quantile(0.25),
        "q3": quantile(0.75),
        "counts": [
            {"level": int(level), "label": level.label, "count": levels.count(level)}
            for level in Level
            if levels.count(level)
        ],
        "description": metric.description,
        "figure": figure_json(fig),
    }
