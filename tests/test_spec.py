from pathlib import Path

import msgspec
import pytest

from qc_report import Level, Spec, SpecError

EXAMPLE = Path(__file__).parents[1] / "examples" / "minimal" / "spec.yaml"


def write(tmp_path, text):
    path = tmp_path / "spec.yaml"
    path.write_text(text)
    return path


def test_example_loads():
    spec = Spec.load(EXAMPLE)
    assert [m.id for m in spec.metrics] == [
        "fragments",
        "aligned_pct",
        "frip",
        "nfr_peak",
    ]
    assert spec.audits[1].severity == Level.ERROR


def test_unknown_field_is_rejected(tmp_path):
    path = write(tmp_path, "title: t\ngroups: []\ncolour: red\n")
    with pytest.raises(msgspec.ValidationError, match="colour"):
        Spec.load(path)


def test_repeated_ids_are_rejected(tmp_path):
    text = (
        "title: t\ngroups:\n  - name: g\n    metrics:\n"
        "      - {id: x, label: X}\n      - {id: x, label: Y}\n"
    )
    with pytest.raises(SpecError, match="repeat"):
        Spec.load(write(tmp_path, text))


def test_audit_on_unknown_metric_is_rejected(tmp_path):
    text = (
        "title: t\ngroups:\n  - name: g\n    metrics:\n      - {id: x, label: X}\n"
        "audits:\n  - {flag: f, level: WARNING, when: y < 1}\n"
    )
    with pytest.raises(SpecError, match="unknown metrics: y"):
        Spec.load(write(tmp_path, text))


def test_unparsable_audit_is_rejected(tmp_path):
    text = (
        "title: t\ngroups:\n  - name: g\n    metrics:\n      - {id: x, label: X}\n"
        "audits:\n  - {flag: f, level: WARNING, when: x <<< 1}\n"
    )
    with pytest.raises(SpecError, match="audit 'f'"):
        Spec.load(write(tmp_path, text))


def test_unknown_level_is_rejected(tmp_path):
    text = (
        "title: t\ngroups:\n  - name: g\n    metrics:\n      - {id: x, label: X}\n"
        "audits:\n  - {flag: f, level: SEVERE, when: x < 1}\n"
    )
    with pytest.raises(msgspec.ValidationError):
        Spec.load(write(tmp_path, text))
