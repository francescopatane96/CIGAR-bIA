import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

def plot_cigar_reads(parsed_reads, chrom, start, end, legend=True, status=None, extra_positions=None):
    """
    Plot deduplicated reads with CIGAR events.
    """

    # Convert DataFrame to list of dicts
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    # Deduplicate
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

    # Draw reads + CIGAR events
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

    # Start / end dashed lines
    ax.axvline(start, color='grey', linestyle='--')
    ax.axvline(end, color='grey', linestyle='--')

    # Extra positions
    if extra_positions:
        for pos in extra_positions[:4]:
            ax.axvline(pos, color='grey', linestyle='--')

    ax.text((start + end) / 2, n_reads - 1, f"{chrom}:{start}-{end}", ha='center', fontsize=9)
    ax.set_xlim(start - 10, end + 10)
    ax.set_ylim(-1, n_reads + 2)
    ax.set_yticks([])
    ax.set_xlabel("Genomic position")
    ax.set_title("Reads with full CIGAR events" + (f" - {status}" if status else ""))

    # Conditional legend
    if legend:
        legend_patches = [
            mpatches.Patch(color='black', label='Match (M)'),
            mpatches.Patch(color='blue', label='Insertion (I)'),
            mpatches.Patch(color='red', label='Deletion (D)'),
            mpatches.Patch(color='purple', label='Skipped region (N)'),
            mpatches.Patch(color='orange', label='Soft clipping (S)')
        ]
        ax.legend(handles=legend_patches, loc='upper left')

    plt.tight_layout()
    plt.show()


import pandas as pd
import matplotlib.pyplot as plt

def plot_editing_distribution(parsed_reads, chrom, start, end, status=None):
    """
    Analizza e mostra distribuzione modifiche e frameshift solo nell'intervallo start-end.
    Considera I, D, N come modifiche per calcolare l'efficienza KO.
    """

    # Converti in lista di dict se necessario
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    # --- Deduplicazione interna e selezione reads che intersecano l'intervallo ---
    seen = set()
    dedup_reads = []
    for read in parsed_reads:
        key = (read.get('barcode'), read.get('umi'))
        if key not in seen and not (read['end'] < start or read['start'] > end):
            seen.add(key)
            dedup_reads.append(read)

    if len(dedup_reads) == 0:
        print("No deduplicated reads in the specified interval.")
        return

    df = pd.DataFrame(dedup_reads)

    # --- Calcolo metriche ---
    total_reads = len(df)
    # considera read editate se I_count, D_count o N_count > 0
    df['edited_any'] = (df['I_count'] + df['D_count'] + df['N_count']) > 0
    edited_reads = df['edited_any'].sum()
    frameshift_reads = df[df['frameshift'] == True].shape[0]

    ko_efficiency = edited_reads / total_reads
    frameshift_fraction = frameshift_reads / total_reads

    print(f"Chromosome {chrom}:{start}-{end}")
    print(f"Total reads: {total_reads}")
    print(f"Edited reads: {edited_reads} -> KO efficiency: {ko_efficiency:.3f}")
    print(f"Reads with frameshift: {frameshift_reads} -> Fraction: {frameshift_fraction:.3f}")

    # --- Prepara dati per stacked barplot ---
    records = []
    for _, row in df.iterrows():
        for edit_type, count in zip(['I','D','N'], [row['I_count'], row['D_count'], row['N_count']]):
            if count > 0:
                records.append({
                    'edit_type': edit_type,
                    'frameshift': row['frameshift'],
                    'reads': 1  # contiamo 1 read per tipo di modifica
                })

    df_long = pd.DataFrame(records)

    if df_long.empty:
        print("No edited reads in the interval.")
        return

    # Pivot per stacked barplot: index = edit_type, colonne = frameshift, valori = numero di reads
    df_plot = df_long.groupby(['edit_type', 'frameshift']).sum().unstack(fill_value=0)['reads']

    # Colori: frameshift True = red, False = green
    colors = {True: 'red', False: 'green'}
    df_plot.plot(kind='bar', stacked=True, color=[colors[col] for col in df_plot.columns], figsize=(6,4))

    plt.ylabel("Number of edited reads")
    plt.xlabel("Edit type")
    plt.title(f"Distribution of edited reads by type and frameshift" + (f" - {status}" if status else ""))
    plt.legend(title="Frameshift")
    plt.tight_layout()
    plt.show()
