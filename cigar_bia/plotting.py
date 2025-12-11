import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

def plot_cigar_reads(parsed_reads, chrom, start, end, status=None, extra_positions=None):
    """
    Plot deduplicated reads with CIGAR events.
    Deduplicazione avviene usando (barcode, UMI).
    
    Args:
        parsed_reads (list of dict OR DataFrame): output di analyze_and_save_editing_events
        chrom (str): Chromosome
        start (int): Start position of target region
        end (int): End position of target region
        status (str, optional): Condition name for plot title
        extra_positions (list of int, optional): up to 4 extra genomic positions da evidenziare con linee tratteggiate
    """

    # Converti DataFrame in lista di dict
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    # Deduplicazione interna
    seen = set()
    dedup_reads = []
    for read in parsed_reads:
        key = (read.get('barcode'), read.get('umi'))
        if key not in seen:
            seen.add(key)
            dedup_reads.append(read)

    n_reads = len(dedup_reads)
    fig_height = min(5, 0.2 * n_reads + 2)
    fig, ax = plt.subplots(figsize=(6, fig_height))

    for y, read in enumerate(dedup_reads):
        ax.hlines(y, read['start'], read['end'], color='black', linewidth=1)
        for op, start_pos, end_or_len in read['cigar']:
            if op == 'I':
                ax.vlines(start_pos, y - 0.3, y + 0.3, color='blue', linewidth=2)
            elif op == 'D':
                ax.hlines(y, start_pos, end_or_len, color='red', linewidth=3)
            elif op == 'N':
                ax.hlines(y, start_pos, end_or_len, color='purple', linewidth=3)
            elif op == 'S':
                ax.vlines(start_pos, y - 0.1, y + 0.1, color='orange', linewidth=2)

    # Linee tratteggiate per start e end
    ax.axvline(start, color='grey', linestyle='--')
    ax.axvline(end, color='grey', linestyle='--')

    # Linee tratteggiate extra se specificate
    if extra_positions:
        for pos in extra_positions[:4]:  # massimo 4 posizioni
            ax.axvline(pos, color='grey', linestyle='--')

    ax.text((start + end) / 2, n_reads - 1, f"{chrom}:{start}-{end}", ha='center', fontsize=9)
    ax.set_xlim(start - 10, end + 10)
    ax.set_ylim(-1, n_reads + 2)
    ax.set_yticks([])
    ax.set_xlabel("Genomic position")
    ax.set_title(f"Reads with full CIGAR events" + (f" - {status}" if status else ""))

    legend_patches = [
        mpatches.Patch(color='black', label='Match (M)'),
        mpatches.Patch(color='blue', label='Insertion (I)'),
        mpatches.Patch(color='red', label='Deletion (D)'),
        mpatches.Patch(color='purple', label='Skipped region (N)'),
        mpatches.Patch(color='orange', label='Soft clipping (S)'),
        mpatches.Patch(color='brown', label='Extra positions')
    ]
    ax.legend(handles=legend_patches, loc='upper left')
    plt.tight_layout()
    plt.show()
