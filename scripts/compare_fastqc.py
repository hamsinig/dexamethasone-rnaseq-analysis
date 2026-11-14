
import argparse
import zipfile
from pathlib import Path
from io import StringIO

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


SAMPLE = "SRR1039508"

MODULES = [
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
]


def extract_fastqc(zip_path):
    """Extract FastQC modules, statuses and tabular data."""

    with zipfile.ZipFile(zip_path, "r") as archive:

        data_file = next(
            name for name in archive.namelist()
            if name.endswith("fastqc_data.txt")
        )

        content = archive.read(data_file).decode("utf-8")

    lines = content.splitlines()
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

            table_text = "\n".join(lines[start:end])

            if table_text.strip():
                df = pd.read_csv(
                    StringIO(table_text),
                    sep="\t",
                    comment=None
                )
            else:
                df = pd.DataFrame()

            modules[module_name] = {
                "status": status.upper(),
                "data": df
            }

            i = end

        i += 1

    return modules


def plot_quality(ax, df):
    """Plot per-base sequence quality."""

    x = range(len(df))

    ax.fill_between(
        x,
        df["Lower Quartile"],
        df["Upper Quartile"],
        alpha=0.2,
        label="25th–75th percentile"
    )

    ax.plot(
        x, df["Mean"],
        label="Mean",
        linewidth=2
    )

    ax.plot(
        x, df["Median"],
        label="Median",
        linestyle="--"
    )

    ax.axhline(
        y=30,
        linestyle=":",
        label="Q30"
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(
        df.iloc[:, 0],
        rotation=45,
        ha="right"
    )

    ax.set_xlabel("Base position")
    ax.set_ylabel("Phred quality score")
    ax.legend(fontsize=8)


def plot_tile_quality(ax, df):
    """Plot per-tile sequence quality as a heatmap."""

    # FastQC uses '#Tile' as the first column name.
    tile_column = (
        "#Tile" if "#Tile" in df.columns else "Tile"
    )

    pivot = df.pivot(
        index=tile_column,
        columns="Base",
        values="Mean"
    )

    image = ax.imshow(
        pivot.values,
        aspect="auto",
        cmap="RdYlBu",
        interpolation="nearest"
    )

    ax.set_ylabel("Tile")
    ax.set_xlabel("Base position")

    plt.colorbar(
        image,
        ax=ax,
        label="Mean quality deviation"
    )


def plot_table(ax, df):
    """Display table-based FastQC modules."""

    ax.axis("off")

    if df.empty:
        ax.text(
            0.5, 0.5,
            "No data available",
            ha="center",
            va="center"
        )
        return

    display_df = df.head(20).fillna("")

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        loc="center",
        cellLoc="center"
    )

    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.3)

    if len(df) > 20:
        ax.set_title(
            f"Showing first 20 of {len(df)} rows"
        )


def plot_generic(ax, df):
    """Plot numeric FastQC data against its first column."""

    if df.empty or len(df.columns) < 2:
        ax.text(
            0.5, 0.5,
            "No graph data available",
            ha="center",
            va="center"
        )
        return

    x_labels = df.iloc[:, 0]
    x = list(range(len(x_labels)))

    plotted = False

    for column in df.columns[1:]:

        values = pd.to_numeric(
            df[column],
            errors="coerce"
        )

        if values.notna().any():
            ax.plot(
                x,
                values,
                label=column,
                linewidth=1.5
            )
            plotted = True

    if not plotted:
        ax.text(
            0.5, 0.5,
            "No numeric data available",
            ha="center",
            va="center"
        )
        return

    ax.set_xticks(x)
    ax.set_xticklabels(
        x_labels,
        rotation=45,
        ha="right"
    )

    ax.set_xlabel(df.columns[0])
    ax.set_ylabel("Value")

    ax.legend(fontsize=7)


def draw_module(ax, module_name, module):
    """Choose the appropriate visualization for each module."""

    if module is None:
        ax.text(
            0.5, 0.5,
            "Module not available",
            ha="center",
            va="center"
        )
        return

    df = module["data"]

    if df.empty:
        ax.text(
            0.5, 0.5,
            "No data available",
            ha="center",
            va="center"
        )
        return

    if module_name in [
        "Basic Statistics",
        "Overrepresented sequences"
    ]:
        plot_table(ax, df)

    elif module_name == "Per base sequence quality":
        plot_quality(ax, df)

    elif module_name == "Per tile sequence quality":
        plot_tile_quality(ax, df)

    else:
        plot_generic(ax, df)

    ax.grid(alpha=0.2)


def create_comparison(raw_dir, trimmed_dir, output_dir):

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    pdf_path = (
        output_dir /
        f"{SAMPLE}_full_fastqc_comparison.pdf"
    )

    with PdfPages(pdf_path) as pdf:

        for read in ["1", "2"]:

            raw_zip = (
                raw_dir /
                f"{SAMPLE}_{read}_fastqc.zip"
            )

            trimmed_zip = (
                trimmed_dir /
                f"{SAMPLE}_{read}_fastqc.zip"
            )

            if not raw_zip.exists():
                raise FileNotFoundError(
                    f"Raw FastQC report not found: {raw_zip}"
                )

            if not trimmed_zip.exists():
                raise FileNotFoundError(
                    f"Trimmed FastQC report not found: {trimmed_zip}"
                )

            print(f"\nProcessing Read {read}...")

            raw_modules = extract_fastqc(raw_zip)
            trimmed_modules = extract_fastqc(trimmed_zip)

            for module_name in MODULES:

                fig, axes = plt.subplots(
                    1, 2,
                    figsize=(16, 7)
                )

                for ax, modules, label in [
                    (axes[0], raw_modules, "RAW"),
                    (axes[1], trimmed_modules, "TRIMMED")
                ]:

                    module = modules.get(module_name)

                    status = (
                        module["status"]
                        if module is not None
                        else "NOT AVAILABLE"
                    )

                    ax.set_title(
                        f"{label} — {status}",
                        fontsize=11
                    )

                    draw_module(
                        ax,
                        module_name,
                        module
                    )

                fig.suptitle(
                    f"{SAMPLE} | Read {read} | {module_name}",
                    fontsize=14
                )

                fig.tight_layout(rect=[0, 0, 1, 0.94])

                pdf.savefig(
                    fig,
                    bbox_inches="tight"
                )

                plt.close(fig)

                print(f"  Completed: {module_name}")

    print("\nComparison PDF successfully created:")
    print(pdf_path)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Compare all FastQC modules "
            "before and after trimming."
        )
    )

    parser.add_argument(
        "--raw",
        required=True,
        help="Directory containing raw FastQC ZIP files"
    )

    parser.add_argument(
        "--trimmed",
        required=True,
        help="Directory containing trimmed FastQC ZIP files"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Directory for comparison PDF"
    )

    args = parser.parse_args()

    create_comparison(
        Path(args.raw),
        Path(args.trimmed),
        Path(args.output)
    )


if __name__ == "__main__":
    main()

