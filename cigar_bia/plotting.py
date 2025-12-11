import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

def plot_cigar_reads(parsed_reads, chrom, start, end, status=None):
    """
    Plot deduplicated reads with CIGAR events.
    La deduplicazione avviene all'interno della funzione usando (barcode, UMI).
    
    Args:
        parsed_reads (list of dict OR DataFrame): output di analyze_and_save_editing_events
        chrom (str): Chromosome
        start (int): Start position of target region
        end (int): End position of target region
        status (str, optional): Condition name for plot title
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

    # Target window
    ax.axvline(start, color='grey', linestyle='--')
    ax.axvline(end, color='grey', linestyle='--')

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
    ]
    ax.legend(handles=legend_patches, loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_editing_summary(parsed_reads, chrom, start, end, status=None):
    """
    Deduplica le reads e mostra grafici di editing:
    1) Efficienza KO (reads editate / tot)
    2) Stacked barplot: distribuzione tipi di editing (I, D, S) e basi modificate
       con colore che indica frameshift o no.
    
    Args:
        parsed_reads (list of dict or DataFrame): output di analyze_and_save_editing_events
        chrom (str): cromosoma
        start (int): inizio regione
        end (int): fine regione
        status (str, optional): nome condizione
    """

    # Converti DataFrame in lista di dict
    if isinstance(parsed_reads, pd.DataFrame):
        parsed_reads = parsed_reads.to_dict(orient='records')

    # Deduplicazione
    seen = set()
    dedup_reads = []
    for read in parsed_reads:
        key = (read.get('barcode'), read.get('umi'))
        if key not in seen:
            seen.add(key)
            dedup_reads.append(read)

    # Trasformiamo in DataFrame per comodità
    df = pd.DataFrame(dedup_reads)

    # --- 1. Efficienza KO ---
    total_reads = len(df)
    edited_reads = df['edited'].sum()
    efficiency = edited_reads / total_reads if total_reads > 0 else 0

    print(f"Chromosome {chrom}:{start}-{end}")
    print(f"Total reads: {total_reads}")
    print(f"Edited reads: {edited_reads}")
    print(f"KO efficiency: {efficiency:.3f}")

    # --- 2. Prepara dati per stacked barplot ---
    # Trasforma in formato lungo
    df_long = df.melt(
        id_vars=['read_name', 'frameshift', 'edited'],
        value_vars=['I_count', 'D_count', 'S_count'],
        var_name='edit_type',
        value_name='count'
    )

    # Filtra solo reads editate
    df_long = df_long[df_long['edited'] == True]

    # Se non ci sono reads editate, fermati
    if df_long.empty:
        print("No edited reads in the specified region.")
        return

    # Stacked barplot: basi modificate per tipo
    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=df_long,
        x='read_name',
        y='count',
        hue='edit_type',
        dodge=False
    )
    plt.xticks([], [])  # nascondi read_name per chiarezza
    plt.ylabel("Number of edited bases")
    plt.xlabel("Reads (deduplicated)")
    plt.title(f"Distribution of editing types per read" + (f" - {status}" if status else ""))
    plt.legend(title="Edit type")
    plt.show()

    # --- 3. Stacked barplot frameshift vs non-frameshift ---
    # Raggruppa totale basi modificate per read
    df_shift = df_long.groupby(['read_name', 'frameshift']).agg({'count': 'sum'}).reset_index()

    # Pivot per stacked barplot
    df_pivot = df_shift.pivot(index='read_name', columns='frameshift', values='count').fillna(0)

    # Colori: frameshift True = red, False = green
    colors = {True: 'red', False: 'green'}
    df_pivot.plot(kind='bar', stacked=True, color=[colors[col] for col in df_pivot.columns], figsize=(8,5))
    plt.ylabel("Number of edited bases")
    plt.xlabel("Reads (deduplicated)")
    plt.title(f"Frameshift vs Non-frameshift" + (f" - {status}" if status else ""))
    plt.legend(title='Frameshift')
    plt.show()

