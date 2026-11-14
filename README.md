# qc-report

One HTML QC report per pipeline run. It grades each sample with ENCODE's audit levels and marks cohort outliers. All figures are interactive plotly figures in a poikilos theme.

## Description

A pipeline hands the engine two things:

- **A metrics table**: a TSV with a `sample_id` column and one column per metric.
- **A spec**: a YAML file that lays out the metrics and lists the audits.

The engine renders one page with these parts:

- A status bar that counts samples at each level: OK, WARNING, NOT COMPLIANT and ERROR. A sample's level is its highest fired audit, as on the ENCODE portal.
- A cohort card per headline metric: the median, the IQR and a histogram stacked by level, with the audit cutoffs drawn in.
- The audit flags that fired, with sample counts. Clicking a flag filters the table.
- A sortable sample table. Each cell is tinted by the highest level among the audits that read it. A blue ▼ or ▲ marks a cohort outlier.
- A detail view per sample: every metric with its level and cohort percentile, the fired flags with their conditions, links, and figures.

Colors follow the portal: green for OK, yellow for WARNING, orange for NOT COMPLIANT, red for ERROR. Every level also prints as text.

### Audits

Each audit is a flag, a level and a condition. The condition is a polars SQL expression over metric ids, so it can copy ENCODE's operators exactly and combine metrics:

```yaml
- flag: moderate number of reproducible peaks
  level: WARNING
  when: overlap_peaks <= 150000 AND idr_peaks <= 70000 AND (overlap_peaks >= 100000 OR idr_peaks >= 50000)
  source: analysis.py:2088-2107
```

Cutoffs MUST come from `~/Projects/encode-data-standards`, the assay's `portal-audits.md`. Copy the flag name, level and operators. Put every deviation in `source`. A null metric never fires an audit.

### Cohort outliers

A metric with `outlier: low`, `high` or `both` flags samples whose modified z-score, `0.6745 * (x - median) / MAD`, passes 3.5 on that side. `scale: log` scores log10 values, which suits counts. Runs under 10 samples get no outliers. Outliers never change a sample's level.

### Figures

Cohort histograms come from the metrics table. Per-sample figures come from the pipeline. The engine ships builders for common tool outputs:
- `read_quality_figure` for fastp
- `gc_bias_figure` for Picard GC bias metrics
- `insert_size_figure` for Picard's insert size histogram, with shaded ranges such as the nucleosome-free one
- `tss_profile_figure` for a TSS enrichment profile: `distance` and `enrichment` columns
- `fingerprint_figure` for a deepTools fingerprint: `fraction_of_bins` and `fraction_of_reads` columns
- `idr_signal_figure` for IDR 2.x output: each pseudoreplicate's signal per peak, split by IDR

`write_sample_panels(path, sample, panels)` writes `samples/<sample>.js` beside the report. A panel is a plotly figure, or the path to a PNG relative to the report. The page loads a sample's file when its detail opens, so a 767-sample run stays light. This works from `file://`.

The page holds plotly.js and the theme's template once. `poikilos.plotly.template(theme)` builds the template from the same theme data as `poikilos.use(theme)`.

## Getting Started

### Dependencies

- Python 3.11 or later.
- polars, plotly, minijinja, msgspec, PyYAML, python-zlib-ng and `poikilos[plotly]`.

### Installing

In a pixi workspace of its own, so report updates never touch the analysis envs:

```toml
[dependencies]
python = "3.14.*"

[pypi-dependencies]
qc-report = { git = "https://github.com/ramirc0/qc-report", tag = "v0.1.0" }
```

qc-report pins poikilos by tag, so pixi installs it too. Pin Python to a
version that every dependency ships wheels for.

### Executing program

```bash
qc-report \
    --spec workflow/report/qc.yaml \
    --metrics results/<run>/report/qc_metrics.tsv \
    --run-dir results/<run> \
    --run-id <run> \
    --version "$(git describe --always --dirty)" \
    --output results/<run>/report/qc_report.html
```

`examples/minimal/` holds a spec and table to start from.

### Adding a pipeline

1. Write a collector script that reads each sample's QC files into one row and writes the metrics table. Use `qc_report.readers` for fastp, cutadapt, samtools flagstat, Picard and deepTools outputs. Parse pipeline-specific files in the script.
2. Write the spec in `workflow/report/qc.yaml`.
   - Group the metrics.
   - Mark the table columns and headline cards.
   - Copy the audits from encode-data-standards.
3. Write a figures script that calls `write_sample_panels` for each sample. Keep the data behind each plot as a pipeline output, so the figure can be drawn from it.
4. Add the collector, figures and `qc-report` rules. Give them their own pixi workspace, so the report's dependencies never change the lock file of the analysis environments.

atac-smk is the reference: `workflow/scripts/gather_qc.py`, `workflow/scripts/qc_figures.py` and `workflow/report/qc.yaml`.

### Spec reference

| Key | Meaning |
| --- | --- |
| `title`, `theme` | Page title, and the poikilos theme for figures and fonts. Default `zed`. |
| `standards` | `name`, `url` and `source` of the audited standard. |
| `groups[].metrics[]` | `id`, `label`, `format` (`decimal`, `count`, `percent`, `flag`), `digits`, `outlier`, `scale`, `table`, `headline`, `description`. |
| `audits[]` | `flag`, `level` (`WARNING`, `NOT_COMPLIANT`, `ERROR`), `when`, `source`. |
| `links[]` | `label` and `path`, relative to the run directory, with `{sample}`. |

Unknown keys, repeated ids and audits that read unknown metrics are errors.

## Help

```bash
uv sync
uv run pytest                      # unit tests and doctests
uv run ruff check && uv run ruff format --check
uv run numpydoc lint src/qc_report/*.py
```

## Authors

Christian Ramirez

## Version History

See the [release history](https://github.com/ramirc0/qc-report/releases).

## Acknowledgments

- Audit levels, colors and cutoffs: [ENCODE-DCC/encoded](https://github.com/ENCODE-DCC/encoded), via encode-data-standards.
- The modified z-score: Iglewicz and Hoaglin, *How to Detect and Handle Outliers*, 1993.
