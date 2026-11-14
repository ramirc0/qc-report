"""Doctests get pytest's tmp_path, so examples can write the files they read."""

import pytest


@pytest.fixture(autouse=True)
def _doctest_tmp_path(doctest_namespace, tmp_path):
    doctest_namespace["tmp_path"] = tmp_path
