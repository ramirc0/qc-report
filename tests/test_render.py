import html as markup
import json
from pathlib import Path

import plotly.graph_objects as go
import polars as pl

from qc_report import Spec, render_report, write_sample_panels
from qc_report.cli import main

EXAMPLE = Path(__file__).parents[1] / "examples" / "minimal"


def test_example_renders(tmp_path):
    output = tmp_path / "report" / "qc_report.html"
    main(
        [
            "--spec", str(EXAMPLE / "spec.yaml"),
            "--metrics", str(EXAMPLE / "metrics.tsv"),
            "--run-dir", str(tmp_path),
            "--run-id", "example",
            "--output", str(output),
        ]
    )  # fmt: skip
    # minijinja escapes "/" too, so compare decoded markup.
    html = markup.unescape(output.read_text())
    assert html.count('<tr data-sample="') == 12
    assert 'id="s-s05"' in html
    # s05 fails depth (ERROR), alignment, FRiP and the NFR peak.
    assert 'data-sample="s05" data-level="3"' in html
    assert "extremely low read depth" in html
    assert 'href="../s05/fastp/s05.html"' in html
    assert "Plotly" in html


def test_sample_names_are_escaped():
    spec = Spec.load(EXAMPLE / "spec.yaml")
    metrics = pl.read_csv(EXAMPLE / "metrics.tsv", separator="\t").with_columns(
        sample_id=pl.lit("<b>x</b>") + pl.col("sample_id")
    )
    html = render_report(metrics, spec, "run")
    assert "<b>x</b>s01" not in html
    assert "&lt;b&gt;x&lt;" in html


def test_sample_panels_carry_figures_and_images(tmp_path):
    path = tmp_path / "samples" / "s1.js"
    fig = go.Figure(go.Scatter(y=[1, 2]), layout={"template": "plotly_dark"})
    write_sample_panels(path, "s1", [("Curve", fig), ("Old plot", "../s1/plot.png")])
    text = path.read_text()
    assert text.startswith("qcReport.addSample(") and text.endswith(");\n")
    payload = json.loads(text.removeprefix("qcReport.addSample(").removesuffix(");\n"))
    assert payload["sample"] == "s1"
    curve, image = payload["panels"]
    assert curve["figure"]["data"][0]["y"] == [1, 2]
    assert "template" not in curve["figure"]["layout"]
    assert image == {"title": "Old plot", "image": "../s1/plot.png"}
