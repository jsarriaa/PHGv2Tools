import argparse
import csv
from typing import List, Dict, Any


HBED_HEADER = [
    "#chrom",
    "start",
    "end",
    "strand",
    "checksum",
    "genome",
    "ref_chr",
    "ref_start",
    "ref_end",
    "ref_checksum",
]


def read_hbed(path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return rows

        for line_no, cols in enumerate(reader, start=2):
            if not cols or len(cols) < 10:
                continue
            try:
                row = {
                    "chrom": cols[0],
                    "start": int(cols[1]),
                    "end": int(cols[2]),
                    "strand": cols[3],
                    "checksum": cols[4],
                    "genome": cols[5],
                    "ref_chr": cols[6],
                    "ref_start": int(cols[7]),
                    "ref_end": int(cols[8]),
                    "ref_checksum": cols[9],
                    "line_no": line_no,
                }
            except ValueError as exc:
                raise ValueError(f"Invalid integer field at line {line_no}: {exc}") from exc

            rows.append(row)

    return rows


def collapse_inversions(rows: List[Dict[str, Any]], max_gap: int = 0) -> List[Dict[str, Any]]:
    inversions = [r for r in rows if r["strand"] == "-"]
    inversions.sort(key=lambda r: (r["chrom"], r["genome"], r["ref_chr"], r["start"], r["end"]))

    collapsed: List[Dict[str, Any]] = []
    current = None

    for row in inversions:
        if current is None:
            current = {
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": "-",
                "genome": row["genome"],
                "ref_chr": row["ref_chr"],
                "ref_start_min": min(row["ref_start"], row["ref_end"]),
                "ref_end_max": max(row["ref_start"], row["ref_end"]),
                "checksums": [row["checksum"]],
                "ref_checksums": [row["ref_checksum"]],
                "blocks": 1,
            }
            continue

        same_group = (
            row["chrom"] == current["chrom"]
            and row["genome"] == current["genome"]
            and row["ref_chr"] == current["ref_chr"]
        )
        concatenated = row["start"] <= (current["end"] + max_gap + 1)

        if same_group and concatenated:
            current["end"] = max(current["end"], row["end"])
            current["ref_start_min"] = min(current["ref_start_min"], row["ref_start"], row["ref_end"])
            current["ref_end_max"] = max(current["ref_end_max"], row["ref_start"], row["ref_end"])
            current["checksums"].append(row["checksum"])
            current["ref_checksums"].append(row["ref_checksum"])
            current["blocks"] += 1
        else:
            collapsed.append(current)
            current = {
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": "-",
                "genome": row["genome"],
                "ref_chr": row["ref_chr"],
                "ref_start_min": min(row["ref_start"], row["ref_end"]),
                "ref_end_max": max(row["ref_start"], row["ref_end"]),
                "checksums": [row["checksum"]],
                "ref_checksums": [row["ref_checksum"]],
                "blocks": 1,
            }

    if current is not None:
        collapsed.append(current)

    normalized: List[Dict[str, Any]] = []
    for idx, group in enumerate(collapsed, start=1):
        normalized.append(
            {
                "chrom": group["chrom"],
                "start": group["start"],
                "end": group["end"],
                "strand": "-",
                "checksum": f"collapsed_{idx}_n{group['blocks']}",
                "genome": group["genome"],
                "ref_chr": group["ref_chr"],
                "ref_start": group["ref_start_min"],
                "ref_end": group["ref_end_max"],
                "ref_checksum": f"collapsed_{idx}_n{group['blocks']}",
                "source_blocks": group["blocks"],
            }
        )

    return normalized


def write_hbed(path: str, rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HBED_HEADER)
        for row in rows:
            writer.writerow(
                [
                    row["chrom"],
                    row["start"],
                    row["end"],
                    row["strand"],
                    row["checksum"],
                    row["genome"],
                    row["ref_chr"],
                    row["ref_start"],
                    row["ref_end"],
                    row["ref_checksum"],
                ]
            )


def summarize(rows: List[Dict[str, Any]], collapsed: List[Dict[str, Any]]) -> List[str]:
    inversion_rows = [r for r in rows if r["strand"] == "-"]
    inversion_lengths = [max(0, r["end"] - r["start"]) for r in inversion_rows]

    by_chrom: Dict[str, Dict[str, int]] = {}
    for row in collapsed:
        chrom = row["chrom"]
        if chrom not in by_chrom:
            by_chrom[chrom] = {"count": 0, "bp": 0}
        by_chrom[chrom]["count"] += 1
        by_chrom[chrom]["bp"] += max(0, row["end"] - row["start"])

    lines = [
        "Inversion collapse summary",
        "==========================",
        f"Input rows: {len(rows)}",
        f"Inversion rows (- strand): {len(inversion_rows)}",
        f"Collapsed inversion rows: {len(collapsed)}",
        f"Total inversion bp (raw): {sum(inversion_lengths):,}",
    ]

    lines.append("")
    lines.append("Per chromosome:")
    for chrom in sorted(by_chrom):
        lines.append(
            f"  {chrom}: collapsed_inversions={by_chrom[chrom]['count']}, collapsed_bp={by_chrom[chrom]['bp']:,}"
        )
    return lines


def write_summary(path: str, lines: List[str]) -> None:
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def plot_collapsed_inversions(rows: List[Dict[str, Any]], collapsed: List[Dict[str, Any]], plot_file: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for --plot-file. Install with: pip install matplotlib"
        ) from exc

    chrom_max: Dict[str, int] = {}
    for row in rows:
        chrom_max[row["chrom"]] = max(chrom_max.get(row["chrom"], 0), row["end"])

    if not chrom_max:
        fig, ax = plt.subplots(figsize=(10, 2))
        ax.text(0.5, 0.5, "No rows found in input h.bed", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(plot_file, dpi=200)
        plt.close(fig)
        return

    by_chrom: Dict[str, List[Dict[str, Any]]] = {}
    for row in collapsed:
        by_chrom.setdefault(row["chrom"], []).append(row)

    chromosomes = sorted(chrom_max.keys())
    fig_height = max(2.0, 0.55 * len(chromosomes) + 1.2)
    fig, ax = plt.subplots(figsize=(13, fig_height))

    global_max = max(chrom_max.values())
    for idx, chrom in enumerate(chromosomes):
        y = len(chromosomes) - idx
        xmax = chrom_max[chrom]

        ax.hlines(y, 0, xmax, color="#b3e5fc", linewidth=6, zorder=1)
        for inv in by_chrom.get(chrom, []):
            ax.hlines(y, inv["start"], inv["end"], color="#d32f2f", linewidth=6, zorder=2)

        ax.text(-0.01 * global_max, y, chrom, ha="right", va="center", fontsize=9)

    ax.set_title("Collapsed inversion blocks (- strand)")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_yticks([])
    ax.set_ylim(0.3, len(chromosomes) + 0.7)
    ax.set_xlim(0, global_max * 1.01)

    for side in ["left", "right", "top"]:
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    fig.savefig(plot_file, dpi=300)
    plt.close(fig)


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Keep only inversion blocks (- strand) from h.bed and collapse concatenated blocks."
    )
    parser.add_argument("input_hbed", help="Input h.bed file")
    parser.add_argument("output_hbed", help="Output collapsed inversion h.bed file")
    parser.add_argument(
        "--max-gap",
        type=int,
        default=0,
        help="Maximum allowed gap (bp) to still merge two inversion blocks (default: 0)",
    )
    parser.add_argument(
        "--summary-file",
        nargs="?",
        const="AUTO",
        default=None,
        help="Optional summary text file path. If used without value, writes <output_hbed>.summary.txt",
    )
    parser.add_argument(
        "--plot-file",
        nargs="?",
        const="AUTO",
        default=None,
        help="Optional output PNG path for inversion haplopaint-style plot. If used without value, writes <output_hbed>.png",
    )

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    if args.max_gap < 0:
        raise ValueError("--max-gap must be >= 0")

    rows = read_hbed(args.input_hbed)
    collapsed = collapse_inversions(rows, max_gap=args.max_gap)
    write_hbed(args.output_hbed, collapsed)

    summary_lines = summarize(rows, collapsed)
    for line in summary_lines:
        print(line)

    print(f"Output written to: {args.output_hbed}")

    if args.summary_file is not None:
        summary_file = (
            f"{args.output_hbed}.summary.txt" if args.summary_file == "AUTO" else args.summary_file
        )
        write_summary(summary_file, summary_lines)
        print(f"Summary written to: {summary_file}")

    if args.plot_file is not None:
        plot_file = f"{args.output_hbed}.png" if args.plot_file == "AUTO" else args.plot_file
        plot_collapsed_inversions(rows, collapsed, plot_file)
        print(f"Plot written to: {plot_file}")


if __name__ == "__main__":
    main()
