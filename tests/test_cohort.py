import polars as pl
import pytest

from qc_report import Group, Metric, Spec, find_outliers, percentiles


def spec(**metric):
    return Spec("QC", [Group("All", [Metric("x", "X", **metric)])])


def frame(values):
    return pl.DataFrame(
        {"sample_id": [f"s{i}" for i in range(len(values))], "x": values}
    )


def test_flags_only_the_bad_side():
    values = [100.0 + i for i in range(18)] + [10.0, 1000.0]
    low = find_outliers(frame(values), spec(outlier="low"))
    assert low.select("sample_id", "side").rows() == [("s18", "low")]
    both = find_outliers(frame(values), spec(outlier="both"))
    assert sorted(both["sample_id"]) == ["s18", "s19"]


def test_log_scale_scores_orders_of_magnitude():
    values = [1e7 * (1 + i / 20) for i in range(19)] + [1e5]
    outliers = find_outliers(frame(values), spec(outlier="low", scale="log"))
    assert outliers["sample_id"].to_list() == ["s19"]


def test_no_spread_flags_nothing():
    assert find_outliers(frame([5.0] * 20), spec(outlier="both")).is_empty()


def test_small_runs_flag_nothing():
    assert find_outliers(frame([1.0, 2.0, 100.0]), spec(outlier="both")).is_empty()


def test_percentiles_rank_within_the_run():
    ranks = percentiles(frame([3.0, 1.0, 2.0, None]), spec())
    assert ranks["x"].to_list()[:3] == pytest.approx([100, 100 / 3, 200 / 3])
    assert ranks["x"][3] is None
