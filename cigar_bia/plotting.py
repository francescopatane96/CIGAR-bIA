import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

def plot_cigar_reads(parsed_reads, chrom, start, end, status=None):
    """
    Plot deduplicated reads con eventi CIGAR.
    Compatibile con l'output di analyze_and_save_editing_events.
    """

    # supporta DataFrame o lista di dict
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    n_reads = len(parsed_reads)
    fig_height = min(5, 0.2 * n_reads + 2)
    fig, ax = plt.subplots(figsize=(6, fig_height))

    for y, read in enumerate(parsed_reads):
        color_line = 'green' if read.get('edited', False) else 'black'
        ax.hlines(y, read['start'], read['end'], color=color_line, linewidth=1)

        for op, start_pos, end_or_len in read.get('cigar', []):
            if op == 'I':
                ax.vlines(start_pos, y - 0.3, y + 0.3, color='blue', linewidth=2)
            elif op == 'D':
                ax.hlines(y, start_pos, end_or_len, color='red', linewidth=3)
            elif op == 'N':
                ax.hlines(y, start_pos, end_or_len, color='purple', linewidth=3)
            elif op == 'S':
                ax.vlines(start_pos, y - 0.1, y + 0.1, color='orange', linewidth=2)

        if read.get('edited', False):
            ax.text(read['end'] + 1, y, f"{read['edited_bases']}bp{' FS' if read['frameshift'] else ''}",
                    fontsize=6, va='center', color='green')

    ax.axvline(start, color='grey', linestyle='--')
    ax.axvline(end, color='grey', linestyle='--')

    ax.text((start + end) / 2, n_reads - 1, f"{chrom}:{start}-{end}", ha='center', fontsize=9)
    ax.set_xlim(start - 10, end + 10)
    ax.set_ylim(-1, n_reads + 2)
    ax.set_yticks([])
    ax.set_xlabel("Genomic position")
    ax.set_title(f"Reads with full CIGAR events" + (f" - {status}" if status else ""))

    legend_patches = [
        mpatches.Patch(color='black', label='Non-edited read'),
        mpatches.Patch(color='green', label='Edited read'),
        mpatches.Patch(color='blue', label='Insertion (I)'),
        mpatches.Patch(color='red', label='Deletion (D)'),
        mpatches.Patch(color='purple', label='Skipped region (N)'),
        mpatches.Patch(color='orange', label='Soft clipping (S)'),
    ]
    ax.legend(handles=legend_patches, loc='upper left')
    plt.tight_layout()
    plt.show()
