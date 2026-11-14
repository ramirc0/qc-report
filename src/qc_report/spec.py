"""The report spec: metric layout, audits, links and theme, read from YAML."""

from enum import IntEnum
from pathlib import Path
from typing import Literal

import msgspec
import polars as pl


class Level(IntEnum):
    """An audit level, in ENCODE's order of severity.

    ``OK`` means no audit fired. The others are the ENCODE portal's levels.

    Examples
    --------
    >>> max(Level.WARNING, Level.NOT_COMPLIANT).name
    'NOT_COMPLIANT'
    """

    OK = 0
    WARNING = 1
    NOT_COMPLIANT = 2
    ERROR = 3

    @property
    def label(self):
        """Return the level as the portal prints it.

        Returns
        -------
        str
            The name with spaces for underscores.

        Examples
        --------
        >>> Level.NOT_COMPLIANT.label
        'NOT COMPLIANT'
        """
        return self.name.replace("_", " ")


class SpecError(ValueError):
    """A spec that repeats metric ids or holds an invalid audit.

    Examples
    --------
    >>> issubclass(SpecError, ValueError)
    True
    """


class Metric(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """One column of the metrics table and how the report shows it.

    Attributes
    ----------
    id : str
        Column name in the metrics table.
    label : str
        Short name for headers, cards and the sample detail.
    format : {"decimal", "count", "percent", "flag"}
        How values print. ``percent`` values already run from 0 to 100.
    digits : int
        Decimal places for ``decimal`` and ``percent``.
    outlier : {"none", "low", "high", "both"}
        The side of the cohort on which a sample is flagged as an outlier.
    scale : {"linear", "log"}
        Scale for outlier detection and the cohort histogram. Counts SHOULD
        use ``log``.
    table : bool
        Show a column in the sample table. The sample detail shows every metric.
    headline : bool
        Show a cohort card with a histogram.
    description : str
        What the metric measures and where it comes from.

    Examples
    --------
    >>> Metric("fragments", "Fragments", format="count", outlier="low").scale
    'linear'
    """

    id: str
    label: str
    format: Literal["decimal", "count", "percent", "flag"] = "decimal"
    digits: int = 2
    outlier: Literal["none", "low", "high", "both"] = "none"
    scale: Literal["linear", "log"] = "linear"
    table: bool = False
    headline: bool = False
    description: str = ""


class Group(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """Metrics shown together under one heading.

    Attributes
    ----------
    name : str
        Heading.
    metrics : list of Metric
        Metrics, in display order.

    Examples
    --------
    >>> Group("Enrichment", [Metric("frip", "FRiP")]).metrics[0].id
    'frip'
    """

    name: str
    metrics: list[Metric]


class Audit(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """A flag that fires on a sample when a condition on its metrics holds.

    Attributes
    ----------
    flag : str
        The flag's name, as the ENCODE portal names it.
    level : {"WARNING", "NOT_COMPLIANT", "ERROR"}
        Severity.
    when : str
        A polars SQL expression over metric ids, such as ``"frip < 0.2"``.
        A null metric never fires it.
    source : str
        Where the cutoff comes from.

    Examples
    --------
    >>> audit = Audit("low FRiP score", "NOT_COMPLIANT", "frip < 0.2")
    >>> audit.severity > Level.WARNING
    True
    """

    flag: str
    level: Literal["WARNING", "NOT_COMPLIANT", "ERROR"]
    when: str
    source: str = ""

    @property
    def severity(self):
        """Return the audit's level.

        Returns
        -------
        Level
            The level as an ordered enum.

        Examples
        --------
        >>> Audit("extremely low read depth", "ERROR", "fragments < 15000000").severity
        <Level.ERROR: 3>
        """
        return Level[self.level]

    def condition(self):
        """Parse the condition.

        Returns
        -------
        polars.Expr
            Boolean expression over metric columns.

        Examples
        --------
        >>> import polars as pl
        >>> low = Audit("low FRiP score", "NOT_COMPLIANT", "frip < 0.2").condition()
        >>> pl.DataFrame({"frip": [0.1, 0.4]}).select(low)["frip"].to_list()
        [True, False]
        """
        return pl.sql_expr(self.when)

    def columns(self):
        """Name the metrics the condition reads.

        Returns
        -------
        list of str
            Metric ids, each once, in order of first appearance.

        Examples
        --------
        >>> Audit("x", "WARNING", "frip >= 0.2 AND frip <= 0.3").columns()
        ['frip']
        """
        return list(dict.fromkeys(self.condition().meta.root_names()))


class Standards(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """The standard the audits copy.

    Attributes
    ----------
    name : str
        Standard, e.g. ``"ENCODE4 ATAC-seq audits"``.
    url : str
        Page that states it.
    source : str
        Exact provenance, such as a file and commit.

    Examples
    --------
    >>> Standards("ENCODE4 ATAC-seq audits").url
    ''
    """

    name: str
    url: str = ""
    source: str = ""


class Link(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """A per-sample file the sample detail links to.

    Attributes
    ----------
    label : str
        Link text.
    path : str
        Path relative to the run directory, with ``{sample}`` for the sample.

    Examples
    --------
    >>> Link("fastp report", "{sample}/fastp/{sample}.html").path.format(sample="s1")
    's1/fastp/s1.html'
    """

    label: str
    path: str


class Spec(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """Everything a pipeline declares about its report.

    Attributes
    ----------
    title : str
        Report title.
    groups : list of Group
        Metrics by heading, in display order.
    audits : list of Audit
        Flags and their conditions.
    links : list of Link
        Per-sample files to link.
    theme : str
        poikilos theme for figures and page fonts.
    standards : Standards or None
        The standard the audits copy.

    Examples
    --------
    >>> spec = Spec("QC", [Group("Depth", [Metric("fragments", "Fragments")])])
    >>> spec.theme
    'zed'
    """

    title: str
    groups: list[Group]
    audits: list[Audit] = []
    links: list[Link] = []
    theme: str = "zed"
    standards: Standards | None = None

    @classmethod
    def load(cls, path):
        r"""Read and check a spec file.

        Parameters
        ----------
        path : str or pathlib.Path
            YAML spec.

        Returns
        -------
        Spec
            The checked spec.

        Raises
        ------
        msgspec.ValidationError
            If a field is missing, unknown or of the wrong type.
        SpecError
            If metric ids repeat, or an audit does not parse or reads an
            unknown metric.

        Examples
        --------
        >>> path = tmp_path / "qc.yaml"
        >>> _ = path.write_text(
        ...     "title: QC\n"
        ...     "groups:\n"
        ...     "  - name: Depth\n"
        ...     "    metrics: [{id: fragments, label: Fragments, format: count}]\n"
        ...     "audits:\n"
        ...     "  - {flag: low read depth, level: WARNING, when: fragments < 25000000}\n"
        ... )
        >>> Spec.load(path).audits[0].columns()
        ['fragments']
        """
        spec = msgspec.yaml.decode(Path(path).read_bytes(), type=cls)
        spec.check()
        return spec

    @property
    def metrics(self):
        """List every metric, in display order.

        Returns
        -------
        list of Metric
            Metrics of all groups.

        Examples
        --------
        >>> groups = [Group("A", [Metric("x", "X")]), Group("B", [Metric("y", "Y")])]
        >>> [metric.id for metric in Spec("QC", groups).metrics]
        ['x', 'y']
        """
        return [metric for group in self.groups for metric in group.metrics]

    def check(self):
        """Check metric ids and audit conditions.

        Raises
        ------
        SpecError
            If metric ids repeat, or an audit does not parse or reads an
            unknown metric.

        Examples
        --------
        >>> spec = Spec(
        ...     "QC",
        ...     [Group("A", [Metric("frip", "FRiP")])],
        ...     [Audit("low FRiP score", "NOT_COMPLIANT", "fripp < 0.2")],
        ... )
        >>> spec.check()
        Traceback (most recent call last):
        ...
        qc_report.spec.SpecError: audit 'low FRiP score' reads unknown metrics: fripp
        """
        ids = [metric.id for metric in self.metrics]
        if repeated := sorted({i for i in ids if ids.count(i) > 1}):
            raise SpecError(f"metric ids repeat: {', '.join(repeated)}")
        for audit in self.audits:
            try:
                columns = audit.columns()
            except pl.exceptions.PolarsError as error:
                raise SpecError(f"audit {audit.flag!r}: {error}") from error
            if unknown := sorted(set(columns) - set(ids)):
                raise SpecError(
                    f"audit {audit.flag!r} reads unknown metrics: {', '.join(unknown)}"
                )
