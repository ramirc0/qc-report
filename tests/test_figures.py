import numpy as np
import pytest

from qc_report import Level, format_count, level_histogram


def test_log_histogram_bins_its_own_counts():
    # Plotly's binning on a log axis read the edges as exponents and the size in
    # data units: about a billion bins, which crashed the browser.
    values = [1.2e7, 2.5e7, 3.1e7, 9.8e7, 4.4e7]
    levels = [Level.WARNING, Level.OK, Level.OK, Level.OK, Level.OK]
    fig = level_histogram(values, levels, bounds=[2e7], log=True)
    assert [trace.type for trace in fig.data] == ["bar", "bar"]
    assert fig.layout.xaxis.type is None
    positions = np.concatenate([trace.x for trace in fig.data])
    assert positions.min() >= np.log10(1.2e7) and positions.max() <= np.log10(9.8e7)
    assert sum(sum(trace.y) for trace in fig.data) == len(values)
    assert "20M" in fig.layout.xaxis.ticktext
    assert fig.layout.shapes[0].x0 == pytest.approx(np.log10(2e7))


def test_linear_histogram_keeps_values():
    fig = level_histogram([0.1, 0.2, 0.25, 0.4], [Level.OK] * 4, bounds=[0.2])
    assert fig.layout.xaxis.ticktext is None
    assert fig.layout.shapes[0].x0 == 0.2


def test_histogram_of_nothing_is_empty():
    assert level_histogram([float("nan")], [Level.OK]).data == ()


@pytest.mark.parametrize(
    ("value", "text"),
    [(0, "0"), (999, "999"), (1000, "1k"), (99950, "100k"), (1.5e9, "1.5B")],
)
def test_format_count(value, text):
    assert format_count(value) == text
