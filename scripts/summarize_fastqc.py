
import argparse
import json
import re
import zipfile
from pathlib import Path
from io import StringIO

import pandas as pd


def extract_modules(zip_path):
    """Extract FastQC module statuses and tabular data."""

    with zipfile.ZipFile(zip_path) as archive:
        data_file = next(
            name for name in archive.namelist()
            if name.endswith("fastqc_data.txt")
        )
        lines = archive.read(data_file).decode("utf-8").splitlines()

    modules = {}
    i = 0

    while i < len(lines):

        if lines[i].startswith(">>") and lines[i] != ">>END_MODULE":

            header = lines[i][2:].split("\t")
            module_name = header[0]
            status = header[1] if len(header) > 1 else "NA"

            start = i + 1
            end = start

            while end < len(lines) and lines[end] != ">>END_MODULE":
                end += 1

            module_lines = lines[start:end]

            # Find the actual tabular header.
            header_index = next(
                (
                    j for j, line in enumerate(module_lines)
                    if line.startswith("#") and "\t" in line
                ),
                None
            )

            df = pd.DataFrame()

            if header_index is not None:
                table_text = "\n".join(module_lines[header_index:])

                df = pd.read_csv(
                    StringIO(table_text),
                    sep="\t"
                )

                # Remove the leading # from column names.
                df.columns = [
                    str(column).lstrip("#")
                    for column in df.columns
                ]

            modules[module_name] = {
                "status": status.upper(),
                "data": df,
                "raw_lines": module_lines
            }

            i = end

        i += 1

    return modules


def numeric_column(df, column):
    """Convert a column to numeric values."""

    if column not in df.columns:
        return pd.Series(dtype=float)

    return pd.to_numeric(df[column], errors="coerce").dropna()


def extract_fastqc_metrics(modules):
    """Calculate comparable numerical metrics from FastQC."""

    results = {}

    # Basic statistics
    basic = modules.get("Basic Statistics", {}).get(
        "data", pd.DataFrame()
    )

    if not basic.empty:
        stats = dict(zip(basic["Measure"], basic["Value"]))

        results["Read length"] = {
            "value": stats.get("Sequence length", "NA"),
            "unit": "bp"
        }

        results["Total sequences"] = {
            "value": stats.get("Total Sequences", "NA"),
            "unit": "reads"
        }

        results["GC content"] = {
            "value": stats.get("%GC", "NA"),
            "unit": "%"
        }

    # Per-base sequence quality
    quality = modules.get(
        "Per base sequence quality", {}
    ).get("data", pd.DataFrame())

    if not quality.empty:
        means = numeric_column(quality, "Mean")
        medians = numeric_column(quality, "Median")

        if not means.empty:
            results["Mean per-base quality"] = {
                "value": means.mean(),
                "unit": "Phred"
            }

            results["Minimum per-base mean quality"] = {
                "value": means.min(),
                "unit": "Phred"
            }

            results["Final-base mean quality"] = {
                "value": means.iloc[-1],
                "unit": "Phred"
            }

            results["Positions with mean Q30 or higher"] = {
                "value": int((means >= 30).sum()),
                "unit": "positions"
            }

        if not medians.empty:
            results["Mean per-base median quality"] = {
                "value": medians.mean(),
                "unit": "Phred"
            }

    # Per-sequence GC distribution
    gc = modules.get(
        "Per sequence GC content", {}
    ).get("data", pd.DataFrame())

    if not gc.empty:
        gc_values = numeric_column(gc, "GC Content")
        counts = numeric_column(gc, "Count")

        if not gc_values.empty and counts.sum() > 0:
            weighted_gc = (gc_values * counts).sum() / counts.sum()

            results["Weighted mean GC"] = {
                "value": weighted_gc,
                "unit": "%"
            }

    # Per-base N content
    n_content = modules.get(
        "Per base N content", {}
    ).get("data", pd.DataFrame())

    if not n_content.empty:
        n_values = numeric_column(n_content, "N-Count")

        if not n_values.empty:
            max_index = n_values.idxmax()

            results["Maximum per-base N content"] = {
                "value": n_values.max(),
                "unit": "%"
            }

            results["Position of maximum N content"] = {
                "value": n_content.loc[max_index, "Base"],
                "unit": "base position"
            }

    # Sequence duplication
    duplication = modules.get(
        "Sequence Duplication Levels", {}
    )

    if duplication:
        for line in duplication["raw_lines"]:
            match = re.match(
                r"#Total Deduplicated Percentage\s+([\d.]+)",
                line
            )

            if match:
                results["Deduplicated percentage"] = {
                    "value": float(match.group(1)),
                    "unit": "%"
                }
                break

    # Adapter content
    adapter = modules.get(
        "Adapter Content", {}
    ).get("data", pd.DataFrame())

    if not adapter.empty:
        for column in adapter.columns:
            if column == "Position":
                continue

            values = numeric_column(adapter, column)

            if not values.empty:
                max_index = values.idxmax()

                results[f"Maximum adapter content: {column}"] = {
                    "value": values.max(),
                    "unit": "%"
                }

                results[f"Position of maximum adapter: {column}"] = {
                    "value": adapter.loc[max_index, "Position"],
                    "unit": "base position"
                }

    return results


def get_module_statuses(modules):
    """Return all FastQC module statuses."""

    return {
        name: module["status"]
        for name, module in modules.items()
    }


def compare_metrics(raw, trimmed, read):
    """Create rows comparing raw and trimmed metrics."""

    rows = []

    all_metrics = sorted(set(raw) | set(trimmed))

    for metric in all_metrics:

        raw_value = raw.get(metric, {}).get("value", "NA")
        trimmed_value = trimmed.get(metric, {}).get("value", "NA")

        unit = (
            raw.get(metric, {}).get("unit")
            or trimmed.get(metric, {}).get("unit")
            or ""
        )

        try:
            delta = float(trimmed_value) - float(raw_value)
        except (ValueError, TypeError):
            delta = "NA"

        rows.append({
            "Read": read,
            "Metric": metric,
            "Raw": raw_value,
            "Trimmed": trimmed_value,
            "Change (Trimmed - Raw)": delta,
            "Unit": unit
        })

    return rows


def flatten_fastp_summary(fastp):
    """Extract important fastp statistics."""

    summary = fastp["summary"]
    before = summary["before_filtering"]
    after = summary["after_filtering"]
    filtering = fastp["filtering_result"]
    adapter = fastp["adapter_cutting"]

    metrics = {
        "Fastp version": summary.get("fastp_version"),
        "Sequencing layout": summary.get("sequencing"),

        "Total reads before filtering": before.get("total_reads"),
        "Total reads after filtering": after.get("total_reads"),

        "Total bases before filtering": before.get("total_bases"),
        "Total bases after filtering": after.get("total_bases"),

        "Q20 rate before (%)": before.get("q20_rate", 0) * 100,
        "Q20 rate after (%)": after.get("q20_rate", 0) * 100,

        "Q30 rate before (%)": before.get("q30_rate", 0) * 100,
        "Q30 rate after (%)": after.get("q30_rate", 0) * 100,

        "Mean read 1 length before (bp)": before.get("read1_mean_length"),
        "Mean read 1 length after (bp)": after.get("read1_mean_length"),

        "Mean read 2 length before (bp)": before.get("read2_mean_length"),
        "Mean read 2 length after (bp)": after.get("read2_mean_length"),

        "GC content before (%)": before.get("gc_content", 0) * 100,
        "GC content after (%)": after.get("gc_content", 0) * 100,

        "Reads passed filtering": filtering.get("passed_filter_reads"),
        "Reads failed filtering": filtering.get("low_quality_reads", 0)
            + filtering.get("too_many_N_reads", 0)
            + filtering.get("too_short_reads", 0)
            + filtering.get("too_long_reads", 0),

        "Reads with adapter trimmed": adapter.get("adapter_trimmed_reads"),
        "Bases trimmed due to adapters": adapter.get("adapter_trimmed_bases"),

        "Duplication rate (%)": fastp.get(
            "duplication", {}
        ).get("rate", 0) * 100,

        "Insert size peak (bp)": fastp.get(
            "insert_size", {}
        ).get("peak", 0)
    }

    total_before = before.get("total_reads", 0)
    total_after = after.get("total_reads", 0)

    metrics["Read retention (%)"] = (
        total_after / total_before * 100
        if total_before else 0
    )

    return metrics


def write_markdown_report(
    output_path,
    sample,
    comparison_df,
    statuses,
    fastp_metrics
):
    """Generate a written QC analysis report."""

    lines = [
        f"# Quality Control Report: {sample}",
        "",
        "## 1. Overview",
        "",
        "This report compares raw and trimmed paired-end RNA-seq "
        "reads using FastQC and fastp.",
        "",
        "FastQC was used to assess sequencing quality before and "
        "after trimming. fastp statistics were used to evaluate "
        "read filtering and adapter removal.",
        "",
        "## 2. Fastp filtering summary",
        "",
        "| Metric | Value |",
        "|---|---:|"
    ]

    for metric, value in fastp_metrics.items():
        if isinstance(value, float):
            value = f"{value:,.3f}"
        elif isinstance(value, int):
            value = f"{value:,}"

        lines.append(f"| {metric} | {value} |")

    lines += [
        "",
        "## 3. FastQC numerical comparison",
        "",
        "| Read | Metric | Raw | Trimmed | Change | Unit |",
        "|---|---|---:|---:|---:|---|"
    ]

    for _, row in comparison_df.iterrows():
        def fmt(value):
            if isinstance(value, (int, float)):
                return f"{value:,.4f}"
            return str(value)

        lines.append(
            f"| {row['Read']} | {row['Metric']} | "
            f"{fmt(row['Raw'])} | {fmt(row['Trimmed'])} | "
            f"{fmt(row['Change (Trimmed - Raw)'])} | "
            f"{row['Unit']} |"
        )

    lines += [
        "",
        "## 4. FastQC module status",
        "",
        "| Read | Module | Raw status | Trimmed status |",
        "|---|---|---|---|"
    ]

    for read, status_data in statuses.items():
        for module, values in status_data.items():
            lines.append(
                f"| {read} | {module} | "
                f"{values['Raw']} | {values['Trimmed']} |"
            )

    lines += [
        "",
        "## 5. Interpretation",
        "",
        "### Read quality",
        "",
        "Compare the per-base quality profiles and Q20/Q30 rates. "
        "Improved average quality after filtering does not necessarily "
        "mean that every retained base was individually improved; "
        "fastp also removes entire read pairs that fail filtering.",
        "",
        "### GC content",
        "",
        "Compare the weighted GC percentage before and after trimming. "
        "A substantial change may indicate that filtering altered "
        "the composition of retained reads.",
        "",
        "### Ambiguous bases",
        "",
        "Compare maximum per-base N content. Lower N content indicates "
        "fewer ambiguous base calls in the retained reads.",
        "",
        "### Adapter content",
        "",
        "Compare the maximum adapter percentage for each adapter type. "
        "A reduction is consistent with successful adapter trimming.",
        "",
        "### Sequence duplication",
        "",
        "High duplication is common in bulk RNA-seq because highly "
        "expressed transcripts generate many identical reads. "
        "FastQC duplication estimates should not automatically be "
        "interpreted as PCR duplication. Reads were not deduplicated "
        "in this workflow.",
        "",
        "## 6. Limitations",
        "",
        "- FastQC module statuses are heuristic quality flags, not "
        "statistical significance tests.",
        "- The FastQC deduplicated percentage is not a direct estimate "
        "of PCR duplication.",
        "- fastp and FastQC calculate duplication using different methods.",
        "- This report describes one sequencing sample and should not "
        "be generalized to all samples until their QC is complete.",
        "- Alignment-level quality metrics will provide additional "
        "evidence about mapping and usable reads.",
        "",
        "## 7. Overall assessment",
        "",
        "The sequencing data should be assessed using the combined "
        "FastQC and fastp results. Final suitability for downstream "
        "RNA-seq analysis should also consider alignment rate, "
        "mapping quality, gene assignment rate and sample-level "
        "consistency.",
        "",
        "## 8. Reproducibility",
        "",
        f"Sample: {sample}",
        "",
        "Tools: FastQC and fastp.",
        "",
        "Numerical values in this report were extracted from the "
        "FastQC ZIP reports and fastp JSON report."
    ]

    output_path.write_text("\n".join(lines) + "\n")


def main():

    parser = argparse.ArgumentParser(
        description="Generate numerical FastQC and fastp QC reports."
    )

    parser.add_argument("--sample", default="SRR1039508")
    parser.add_argument("--raw", required=True)
    parser.add_argument("--trimmed", required=True)
    parser.add_argument("--fastp", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    sample = args.sample

    raw_dir = Path(args.raw)
    trimmed_dir = Path(args.trimmed)
    fastp_path = Path(args.fastp)
    output_dir = Path(args.output)

    output_dir.mkdir(parents=True, exist_ok=True)

    comparison_rows = []
    statuses = {}

    for read in ["1", "2"]:

        raw_zip = raw_dir / f"{sample}_{read}_fastqc.zip"
        trimmed_zip = trimmed_dir / f"{sample}_{read}_fastqc.zip"

        print(f"Processing Read {read}...")

        raw_modules = extract_modules(raw_zip)
        trimmed_modules = extract_modules(trimmed_zip)

        raw_metrics = extract_fastqc_metrics(raw_modules)
        trimmed_metrics = extract_fastqc_metrics(trimmed_modules)

        comparison_rows.extend(
            compare_metrics(raw_metrics, trimmed_metrics, f"R{read}")
        )

        raw_status = get_module_statuses(raw_modules)
        trimmed_status = get_module_statuses(trimmed_modules)

        statuses[f"R{read}"] = {}

        for module in sorted(set(raw_status) | set(trimmed_status)):
            statuses[f"R{read}"][module] = {
                "Raw": raw_status.get(module, "NOT AVAILABLE"),
                "Trimmed": trimmed_status.get(module, "NOT AVAILABLE")
            }

    comparison_df = pd.DataFrame(comparison_rows)

    comparison_csv = output_dir / f"{sample}_qc_comparison.csv"
    comparison_df.to_csv(comparison_csv, index=False)

    with open(fastp_path) as f:
        fastp = json.load(f)

    fastp_metrics = flatten_fastp_summary(fastp)

    fastp_df = pd.DataFrame(
        list(fastp_metrics.items()),
        columns=["Metric", "Value"]
    )

    fastp_csv = output_dir / f"{sample}_fastp_summary.csv"
    fastp_df.to_csv(fastp_csv, index=False)

    report_path = output_dir / f"{sample}_qc_report.md"

    write_markdown_report(
        report_path,
        sample,
        comparison_df,
        statuses,
        fastp_metrics
    )

    print("\nQC analysis completed!")
    print(f"Numerical comparison: {comparison_csv}")
    print(f"Fastp summary:        {fastp_csv}")
    print(f"Written report:       {report_path}")


if __name__ == "__main__":
    main()

