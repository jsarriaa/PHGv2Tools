#!/usr/bin/env python3
"""
UpSet plot generator for PHG pangenome hapIDranges.tsv files.

Provides a `main(argv)` entry so it can be called from the `phgtools` CLI.
"""
import sys
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate an UpSet plot from a Pangenome Haplotype Hash table."
    )
    parser.add_argument(
        "input_file",
        help="Path to the input TSV file (e.g., hapIDranges.tsv)"
    )
    parser.add_argument(
        "--mode",
        choices=['blocks', 'bp'],
        default='blocks',
        help="Metric to calculate intersections: 'blocks' (count) or 'bp' (base pairs)."
    )
    parser.add_argument(
        "--top",
        type=int,
        default=30,
        help="Number of top intersections to display (default: 30)."
    )
    parser.add_argument(
        "--out",
        default="pangenome_upset.png",
        help="Output image filename (default: pangenome_upset.png)."
    )
    return parser.parse_args(argv)


def process_data(file_path, mode):
    # Read file, tolerate occasional header typos where two names are squashed
    with open(file_path, 'r') as f:
        lines = f.readlines()

    if len(lines) == 0:
        sys.exit("Error: File is empty.")

    if "HOR_3365HOR_3474" in lines[0]:
        lines[0] = lines[0].replace("HOR_3365HOR_3474", "HOR_3365    HOR_3474")

    import io
    df = pd.read_csv(io.StringIO("".join(lines)), sep=r'\s+', index_col=False)

    metadata_cols = ['#CHROM', 'START', 'END']
    genome_cols = [c for c in df.columns if c not in metadata_cols]

    df['LENGTH'] = df['END'] - df['START']

    combinations = []
    for _, row in df.iterrows():
        weight = 1 if mode == 'blocks' else row['LENGTH']
        valid_cells = row[genome_cols].replace('.', np.nan).dropna()
        unique_hashes = valid_cells.unique()
        for h in unique_hashes:
            members = tuple(sorted(valid_cells[valid_cells == h].index.tolist()))
            combinations.append({'members': members, 'weight': weight})

    comb_df = pd.DataFrame(combinations)
    intersection_groups = comb_df.groupby('members')['weight'].sum().reset_index(name='count')
    intersection_groups = intersection_groups.sort_values(by='count', ascending=False).reset_index(drop=True)

    set_sizes = {g: 0 for g in genome_cols}
    for _, row in intersection_groups.iterrows():
        for g in row['members']:
            set_sizes[g] += row['count']

    return intersection_groups, set_sizes, genome_cols


def plot_upset(intersection_groups, set_sizes, genome_cols, top_n, output_img, mode):
    top_intersections = intersection_groups.head(top_n)
    sorted_genomes = sorted(genome_cols, key=lambda x: set_sizes[x], reverse=False)
    n_sets = len(sorted_genomes)

    counts = top_intersections['count'].tolist()
    matrix = np.zeros((n_sets, len(top_intersections)))

    for j, row in top_intersections.iterrows():
        for i, g in enumerate(sorted_genomes):
            if g in row['members']:
                matrix[i, j] = 1

    fig, axs = plt.subplots(2, 2, figsize=(16, 10), 
                            gridspec_kw={'width_ratios': [1, 4], 'height_ratios': [3, 2],
                                         'wspace': 0.12, 'hspace': 0.05})
    axs[0, 0].axis('off') 
    ax_hist, ax_set, ax_matrix = axs[0, 1], axs[1, 0], axs[1, 1]

    metric_label = "Base Pairs (bp)" if mode == 'bp' else "Number of Blocks"
    y_formatter = lambda x: f"{int(x):,}"

    # Intersection sizes
    ax_hist.bar(range(len(counts)), counts, color='#2c3e50', width=0.6, edgecolor='black', alpha=0.9)
    ax_hist.set_ylabel(f'Intersection Size\n({metric_label})', fontsize=12, fontweight='bold')
    ax_hist.set_xticks([])
    ax_hist.set_xlim(-0.6, len(counts) - 0.2)
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)

    max_val = max(counts) if counts else 1
    ax_hist.set_ylim(0, max_val * 1.35) 
    for i, val in enumerate(counts):
        y_offset = max_val * 0.02
        ax_hist.text(i, val + y_offset, y_formatter(val), 
                     ha='left', va='bottom', fontsize=9, rotation=45)

    # Total genome sizes
    sorted_sizes = [set_sizes[g] for g in sorted_genomes]
    ax_set.barh(range(n_sets), sorted_sizes, color='#34495e', height=0.5, edgecolor='black', alpha=0.9)
    ax_set.set_xlabel(f'Total Size ({metric_label})', fontsize=12, fontweight='bold')
    ax_set.set_yticks(range(n_sets))
    ax_set.set_yticklabels(sorted_genomes, fontsize=11, fontweight='bold')

    min_size = min(sorted_sizes) if sorted_sizes else 0
    max_size = max(sorted_sizes) if sorted_sizes else 1
    size_range = max_size - min_size
    if size_range == 0:
        size_range = max_size * 0.1 if max_size > 0 else 1
    lower_lim = max(0, min_size - (size_range * 0.5))
    upper_lim = max_size + (size_range * 0.15)
    ax_set.set_xlim(upper_lim, lower_lim)
    ax_set.spines['top'].set_visible(False)
    ax_set.spines['left'].set_visible(False)
    ax_set.xaxis.grid(True, linestyle='--', alpha=0.5)
    ax_set.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    plt.setp(ax_set.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Matrix
    ax_matrix.set_xlim(-0.5, len(counts) - 0.5)
    ax_matrix.set_ylim(-0.5, n_sets - 0.5)
    ax_matrix.set_xticks(range(len(counts)))
    ax_matrix.set_xticklabels([])
    ax_matrix.set_yticks(range(n_sets))
    ax_matrix.set_yticklabels([])

    for s in range(n_sets):
        bg_color = '#f8f9fa' if s % 2 == 0 else '#ffffff'
        ax_matrix.axhspan(s - 0.5, s + 0.5, color=bg_color, zorder=0)
        ax_matrix.axhline(s, color='#bdc3c7', linestyle='-', linewidth=0.5, alpha=0.4, zorder=1)

    for j in range(len(counts)):
        active_idx = np.where(matrix[:, j] == 1)[0]
        for s in range(n_sets):
            if matrix[s, j] == 0:
                ax_matrix.plot(j, s, marker='o', markersize=10, color='#dcdde1', zorder=2)
        if len(active_idx) > 1:
            ax_matrix.plot([j, j], [min(active_idx), max(active_idx)], color='#27ae60', linestyle='-', linewidth=3, zorder=3)
        for s in active_idx:
            ax_matrix.plot(j, s, marker='o', markersize=12, color='#27ae60', zorder=4)

    for spine in ['top', 'right', 'left', 'bottom']:
        ax_matrix.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    plt.close()


def main(argv=None):
    args = parse_arguments(argv)
    if not os.path.exists(args.input_file):
        sys.exit(f"Error: File '{args.input_file}' not found.")
    intersections, sizes, genomes = process_data(args.input_file, args.mode)
    plot_upset(intersections, sizes, genomes, args.top, args.out, args.mode)


if __name__ == '__main__':
    main()
