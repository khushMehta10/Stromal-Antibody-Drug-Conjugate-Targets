#!/usr/bin/env python3
"""
Gene Frequency Analysis Across Cancer Subtypes
Analyzes gene occurrence across patients within each subtype
Identifies universal, common, and subtype-specific genes
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('default')
sns.set_palette("husl")

################################################################################
# PARAMETERS
################################################################################

# Known stromal markers to note in plots
known_markers = ['FAP', 'PDPN', 'COL1A1', 'COL1A2', 'COL11A1', 'POSTN', 'THY1',
                 'PDGFRA', 'PDGFRB', 'VIM', 'DCN', 'LUM', 'ACTA2', 'FN1',
                 'S100A4', 'SPARC', 'TAGLN', 'MMP2', 'MMP11', 'VCAN', 'CDH11']

# Frequency thresholds to analyze
FREQUENCY_THRESHOLDS = {
    'Universal': 1.0,      # Present in 100% of patients
    'Frequent': 0.75,   # Present in ≥75% of patients
    'Very_Common': 0.50,        # Present in ≥50% of patients
    'Common': 0.25       # Present in ≥25% of patients
}

# Top N genes to show in figures
TOP_N_GENES = 25

# Minimum number of patients required for analysis
MIN_PATIENTS = 3

# Output directory
OUTPUT_DIR = "gene_frequency_analysis"

################################################################################
# FUNCTIONS
################################################################################

def find_subtype_folders(base_dir):
    """Find all subdirectories (assumed to be subtypes)"""
    base_path = Path(base_dir)
    subtype_dirs = [d for d in base_path.iterdir() if d.is_dir()]
    return sorted(subtype_dirs)

def load_patient_data(subtype_dir):
    """Load all CSV files from a subtype directory"""
    csv_files = list(Path(subtype_dir).glob("*.csv"))
    
    if not csv_files:
        return None, [], None
    
    patient_data = []
    patient_ids = []
    
    # Try to load first file to detect format
    first_file = csv_files[0]
    df_test = None
    separator = '\t'
    
    # Try tab-separated first
    try:
        df_test = pd.read_csv(first_file, sep='\t', nrows=5)
    except:
        pass
    
    # Try comma-separated if tab failed
    if df_test is None or len(df_test.columns) < 3:
        try:
            df_test = pd.read_csv(first_file, sep=',', nrows=5)
            separator = ','
        except:
            pass
    
    if df_test is None:
        print(f"  Error: Could not detect file format for {first_file.name}")
        return None, [], None
    
    # Detect column names
    columns = df_test.columns.tolist()
    print(f"  Detected columns: {columns[:5]}..." if len(columns) > 5 else f"  Detected columns: {columns}")
    
    # Find gene column
    gene_col = None
    for possible_name in ['gene', 'Gene', 'GENE', 'gene_name', 'symbol', 'Symbol']:
        if possible_name in columns:
            gene_col = possible_name
            break
    
    if gene_col is None:
        print(f"  ERROR: No gene column found. Available columns: {columns}")
        print(f"  Please specify which column contains gene names")
        return None, [], None
    
    print(f"  Using '{gene_col}' as gene column")
    
    # Load all files
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file, sep=separator)
            
            # Verify gene column exists
            if gene_col not in df.columns:
                print(f"  Warning: '{gene_col}' column not found in {csv_file.name}, skipping")
                continue
            
            patient_id = csv_file.stem
            patient_data.append(df)
            patient_ids.append(patient_id)
        except Exception as e:
            print(f"  Warning: Could not read {csv_file.name}: {e}")
    
    return patient_data, patient_ids, gene_col

def analyze_gene_frequency(patient_data, patient_ids, gene_col='gene'):
    """
    Analyze gene frequency across patients
    Returns DataFrame with gene statistics
    """
    n_patients = len(patient_data)
    
    # Collect all genes and their metrics
    gene_occurrences = defaultdict(list)
    gene_ranks = defaultdict(list)
    gene_log2fc = defaultdict(list)
    gene_stromal_exclusivity = defaultdict(list)
    
    for patient_df in patient_data:
        genes_in_patient = set(patient_df[gene_col].values)
        
        for gene in genes_in_patient:
            gene_occurrences[gene].append(1)
            
            # Get metrics for this gene
            gene_data = patient_df[patient_df[gene_col] == gene].iloc[0]
            
            # Safely get rank if available
            if 'rank' in patient_df.columns:
                gene_ranks[gene].append(gene_data['rank'])
            
            # Safely get log2fc if available
            if 'log2fc_vs_all' in patient_df.columns:
                gene_log2fc[gene].append(gene_data['log2fc_vs_all'])
            elif 'log2FC' in patient_df.columns:
                gene_log2fc[gene].append(gene_data['log2FC'])
            
            # Safely get stromal exclusivity if available
            if 'stromal_exclusivity' in patient_df.columns:
                gene_stromal_exclusivity[gene].append(gene_data['stromal_exclusivity'])
    
    # Create summary dataframe
    gene_summary = []
    
    for gene in gene_occurrences.keys():
        n_occurrences = len(gene_occurrences[gene])
        frequency = n_occurrences / n_patients
        
        summary_entry = {
            'gene': gene,
            'n_patients': n_occurrences,
            'frequency': frequency,
            'frequency_pct': frequency * 100
        }
        
        # Add optional metrics
        if gene_ranks[gene]:
            summary_entry['mean_rank'] = np.mean(gene_ranks[gene])
            summary_entry['median_rank'] = np.median(gene_ranks[gene])
        
        if gene_log2fc[gene]:
            summary_entry['mean_log2fc'] = np.mean(gene_log2fc[gene])
            summary_entry['median_log2fc'] = np.median(gene_log2fc[gene])
        
        if gene_stromal_exclusivity[gene]:
            summary_entry['mean_stromal_exclusivity'] = np.mean(gene_stromal_exclusivity[gene])
            summary_entry['median_stromal_exclusivity'] = np.median(gene_stromal_exclusivity[gene])
        
        gene_summary.append(summary_entry)
    
    df_summary = pd.DataFrame(gene_summary)
    df_summary = df_summary.sort_values('frequency', ascending=False)

    # Add composite score (higher is better) if rank is available
    if 'mean_rank' in df_summary.columns:
        # mean_rank: lower is better; convert to a higher-better score
        df_summary['composite_score'] = df_summary['frequency_pct'] / (df_summary['mean_rank'].replace(0, np.nan))
        df_summary['composite_score'] = df_summary['composite_score'].fillna(0)
    else:
        df_summary['composite_score'] = df_summary['frequency_pct']

    return df_summary

def categorize_genes(df_summary, n_patients):
    """Categorize genes by frequency thresholds"""
    categorized = {}
    
    for category, threshold in FREQUENCY_THRESHOLDS.items():
        genes = df_summary[df_summary['frequency'] >= threshold]['gene'].tolist()
        categorized[category] = genes
        print(f"  {category} (≥{threshold*100:.0f}%): {len(genes)} genes")
    
    return categorized

def create_frequency_barplot(df_summary, subtype_name, output_dir):
    """Create bar plot of gene frequencies"""
    
    # Top genes by frequency
    top_genes = df_summary.head(TOP_N_GENES).copy()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create barplot
    bars = ax.barh(range(len(top_genes)), top_genes['frequency_pct'], 
                   color="#325ea8")
    
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes['gene'], fontsize=10)
    for lab in  ax.get_yticklabels():
        t = lab.get_text()
        if t in known_markers:
            lab.set_color("crimson")
            lab.set_fontweight("bold")
    ax.set_xlabel('Frequency (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=12, fontweight='bold')
    ax.set_title(f'{subtype_name}: Top {TOP_N_GENES} Genes by Frequency\n'
                 f'(n={top_genes["n_patients"].iloc[0]} patients)',
                 fontsize=14, fontweight='bold')
    
    # Add frequency threshold lines
    for category, threshold in FREQUENCY_THRESHOLDS.items():
        ax.axvline(threshold * 100, color='red', linestyle='--', 
                  alpha=0.3, linewidth=1)
        ax.text(threshold * 100, len(top_genes)-1, f'  {category}', 
               fontsize=8, color='red', rotation=90, va='top')
    
    ax.set_xlim(0, 105)
    ax.grid(axis='x', alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{subtype_name}_frequency_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved frequency barplot")

def _savefig_png(path, dpi=300):
    """Save figure as PNG (consistent output format)."""
    plt.savefig(path, dpi=dpi, bbox_inches='tight')


def create_ranked_barplot(df_summary, subtype_name, output_dir, sort_by, ascending, filename_suffix,
                          xlabel, title, top_n=None):
    """Generic horizontal barplot for ranked gene lists."""
    if top_n is None:
        top_n = TOP_N_GENES

    if sort_by not in df_summary.columns:
        print(f"  Skipping {filename_suffix} (missing column: {sort_by})")
        return

    plot_df = df_summary.sort_values(sort_by, ascending=ascending).head(top_n).copy()

    fig, ax = plt.subplots(figsize=(8, 6))

    # Use viridis mapped to the plotted values (normalized)
    vals = plot_df[sort_by].astype(float).values
    vmin, vmax = np.nanmin(vals), np.nanmax(vals)
    denom = (vmax - vmin) if vmax != vmin else 1.0
    # colors = plt.cm.viridis((vals - vmin) / denom)

    ax.barh(range(len(plot_df)), vals, color='#325ea8', zorder=3)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['gene'], fontsize=10)

    # highlight known markers
    km = set(globals().get('known_markers', []))
    for lab in ax.get_yticklabels():
        if lab.get_text() in km:
            lab.set_color('crimson')
            lab.set_fontweight('bold')

    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=12, fontweight='bold')
    # n = number of samples in the group (a gene's sample count divided by its frequency)
    ax.set_title(title + '\n'
                 f'(n={int(round(100 * plot_df["n_patients"].iloc[0] / plot_df["frequency_pct"].iloc[0]))} samples)',
                 fontsize=14, fontweight='bold')
    
    ax.grid(alpha=0.3)
    ax.grid(True, zorder=0)
    ax.invert_yaxis()  # top item at top

    plt.tight_layout()
    out_path = f"{output_dir}/{subtype_name}_{filename_suffix}.png"
    _savefig_png(out_path)
    plt.close()
    print(f"  Saved {filename_suffix} barplot")


def create_subtype_ranked_plots(df_summary, subtype_name, output_dir):
    """Create separate ranked plots by frequency, mean_rank, and composite score."""
    # Frequency (descending)
    create_ranked_barplot(
        df_summary, subtype_name, output_dir,
        sort_by='frequency_pct', ascending=False,
        filename_suffix='frequency',
        xlabel='Frequency (%)',
        title=f"{subtype_name}: Top {TOP_N_GENES} Genes by Frequency"
    )

    # Mean rank (ascending = better)
    if 'mean_rank' in df_summary.columns:
        create_ranked_barplot(
            df_summary, subtype_name, output_dir,
            sort_by='mean_rank', ascending=True,
            filename_suffix='mean_rank',
            xlabel='Mean Rank (lower is better)',
            title=f"{subtype_name}: Top {TOP_N_GENES} Genes by Mean Rank"
        )

    # Composite (descending)
    if 'composite_score' in df_summary.columns:
        create_ranked_barplot(
            df_summary, subtype_name, output_dir,
            sort_by='composite_score', ascending=False,
            filename_suffix='composite',
            xlabel='Composite Score (freq% / mean_rank)',
            title=f"{subtype_name}: Top {TOP_N_GENES} Genes by Composite Score"
        )


def create_all_subtypes_top_targets_plot(all_subtypes_data, output_dir, top_n=20, rank_metric='composite_score'):
    """Summary plot: same y-axis genes across subtypes, showing top targets per subtype."""
    if not all_subtypes_data:
        return

    # Build a global ranking from subtype summaries
    subtype_summaries = {k: v['summary'].copy() for k, v in all_subtypes_data.items() if 'summary' in v}
    if len(subtype_summaries) < 2:
        return

    # Ensure rank_metric exists
    for st, df in subtype_summaries.items():
        if rank_metric not in df.columns:
            # fall back
            rank_metric = 'frequency_pct'
            break

    # Aggregate: mean of rank_metric across subtypes (missing -> 0)
    all_genes = sorted(set().union(*[set(df['gene']) for df in subtype_summaries.values()]))
    agg = pd.DataFrame({'gene': all_genes})
    for st, df in subtype_summaries.items():
        agg = agg.merge(df[['gene', rank_metric]].rename(columns={rank_metric: st}), on='gene', how='left')
    metric_cols = list(subtype_summaries.keys())
    agg[metric_cols] = agg[metric_cols].fillna(0)

    agg['mean_metric'] = agg[metric_cols].mean(axis=1)
    top_global = agg.sort_values('mean_metric', ascending=False).head(top_n)['gene'].tolist()

    # Plot
    subtypes = list(subtype_summaries.keys())
    n = len(subtypes)
    fig_w = max(10, 4 * n)
    fig_h = max(6, 0.35 * top_n + 2)
    fig, axes = plt.subplots(1, n, figsize=(fig_w, fig_h), sharey=True)

    if n == 1:
        axes = [axes]

    y = np.arange(len(top_global))

    for ax, st in zip(axes, subtypes):
        df = subtype_summaries[st].set_index('gene')
        vals = [float(df.loc[g, rank_metric]) if g in df.index else 0.0 for g in top_global]
        ax.barh(y, vals)
        ax.set_title(f"{st} Top {top_n}", fontsize=12, fontweight='bold')
        ax.set_xlabel(rank_metric)
        ax.grid(alpha=0.25)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(top_global, fontsize=9)
    axes[0].invert_yaxis()

    fig.suptitle(f"Top {top_n} Stromal Targets by Subtype ({rank_metric})", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = f"{output_dir}/All_Subtypes_Top{top_n}_{rank_metric}.png"
    _savefig_png(out_path)
    plt.close()
    print(f"Saved all-subtypes top targets plot → {out_path}")


def create_all_subtypes_top_targets_plot_per_subtype(all_subtypes_data, output_dir, top_n=10, rank_metric='composite_score'):
    """Summary plot: each subtype subplot uses its *own* top-N genes (not shared y-axis).

    Generates a single figure with one subplot per subtype.
    - For frequency/composite: higher is better (descending)
    - For mean_rank: lower is better (ascending)
    """
    if not all_subtypes_data:
        return

    subtype_summaries = {k: v['summary'].copy() for k, v in all_subtypes_data.items() if 'summary' in v}
    if not subtype_summaries:
        return

    # Determine sort direction
    ascending = True if rank_metric == 'mean_rank' else False

    # Validate metric availability; fall back for missing
    for st, df in subtype_summaries.items():
        if rank_metric not in df.columns:
            if rank_metric == 'mean_rank' and 'mean_rank' not in df.columns:
                print(f"  Skipping per-subtype top plot for mean_rank (missing mean_rank in {st})")
                return
            if rank_metric == 'composite_score' and 'composite_score' not in df.columns:
                # composite needs rank; fall back to frequency
                rank_metric = 'frequency_pct'
                ascending = False
                break
            if rank_metric not in df.columns:
                rank_metric = 'frequency_pct'
                ascending = False
                break

    subtypes = list(subtype_summaries.keys())
    n = len(subtypes)
    fig_w = max(10, 4 * n)
    fig_h = max(6, 0.35 * top_n + 2)
    fig, axes = plt.subplots(1, n, figsize=(8, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, st in zip(axes, subtypes):
        df = subtype_summaries[st]
        plot_df = df.sort_values(rank_metric, ascending=ascending).head(top_n).copy()
        # Ensure "best" at top visually
        y = np.arange(len(plot_df))
        ax.barh(y, plot_df[rank_metric].astype(float).values, color='#325ea8', zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(plot_df['gene'].tolist(), fontsize=12)
        ax.invert_yaxis()
        ax.set_title(f"{st} Top {top_n}", fontsize=12, fontweight='bold')
        ax.set_xlabel(rank_metric if rank_metric != 'mean_rank' else 'Mean Rank (lower is better)')
        ax.grid(alpha=0.25, zorder=0)

        # highlight known markers
        km = set(globals().get('known_markers', []))
        for lab in ax.get_yticklabels():
            if lab.get_text() in km:
                lab.set_color('crimson')
                lab.set_fontweight('bold')

    fig.suptitle(f"Per-Subtype Top {top_n} Stromal Targets ({rank_metric})", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = f"{output_dir}/All_Subtypes_PerSubtypeTop{top_n}_{rank_metric}.png"
    _savefig_png(out_path)
    plt.close()
    print(f"Saved per-subtype top targets plot → {out_path}")

def create_rank_vs_frequency_scatter(df_summary, subtype_name, output_dir):
    """Create scatter plot of rank vs frequency"""
    
    # Check if required columns exist
    if 'mean_rank' not in df_summary.columns:
        print(f"  Skipping rank vs frequency plot (no rank data)")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Use log2fc for color if available, otherwise use constant color
    if 'mean_log2fc' in df_summary.columns:
        scatter = ax.scatter(df_summary['mean_rank'], 
                            df_summary['frequency_pct'],
                            c=df_summary['mean_log2fc'],
                            s=100,
                            alpha=0.6,
                            cmap='RdYlBu_r',
                            edgecolors='black',
                            linewidth=0.5)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Mean log2FC vs All', fontsize=11, fontweight='bold')
    else:
        ax.scatter(df_summary['mean_rank'], 
                  df_summary['frequency_pct'],
                  s=100,
                  alpha=0.6,
                  color='steelblue',
                  edgecolors='black',
                  linewidth=0.5)
    
    # Label top genes
    top_genes = df_summary.head(20)
    for idx, row in top_genes.iterrows():
        ax.annotate(row['gene'], 
                   (row['mean_rank'], row['frequency_pct']),
                   fontsize=7,
                   alpha=0.7,
                   xytext=(5, 5),
                   textcoords='offset points')
    
    ax.set_xlabel('Mean Rank Across Patients', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'{subtype_name}: Gene Rank vs Frequency\n'
                 f'Color = Mean log2FC',
                 fontsize=14, fontweight='bold')
    
    # Add frequency threshold lines
    for category, threshold in FREQUENCY_THRESHOLDS.items():
        ax.axhline(threshold * 100, color='gray', linestyle='--', 
                  alpha=0.3, linewidth=1)
    
    ax.grid(alpha=0.3)

    # Most frequent at the top
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{output_dir}/{subtype_name}_rank_vs_frequency.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved rank vs frequency plot")

def create_presence_absence_heatmap(patient_data, patient_ids, subtype_name,
                                    output_dir, gene_col='gene', top_n=25):
    """Create heatmap showing gene presence/absence across patients.
       Renames patient IDs starting with 'CID' to 'FID' for sorting/labels.
    """

    from collections import defaultdict
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    # -------------------------------------------------
    # 1) Rename CID -> FID (for display + sorting)
    # -------------------------------------------------
    renamed_patient_ids = [
        ("FID" + pid[3:]) if isinstance(pid, str) and pid.startswith("CID") else pid
        for pid in patient_ids
    ]

    # Sort patients using the renamed IDs, but keep data aligned
    sorted_pairs = sorted(zip(renamed_patient_ids, patient_ids, patient_data),
                          key=lambda x: x[0])
    sorted_renamed_ids, sorted_original_ids, sorted_patient_data = zip(*sorted_pairs)

    # -------------------------------------------------
    # 2) Get gene frequencies
    # -------------------------------------------------
    gene_counts = defaultdict(int)
    for patient_df in sorted_patient_data:
        for gene in patient_df[gene_col].values:
            gene_counts[gene] += 1

    # Select top genes
    top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_gene_names = [g[0] for g in top_genes]

    # -------------------------------------------------
    # 3) Build presence/absence matrix
    # -------------------------------------------------
    matrix = []
    for patient_df in sorted_patient_data:
        patient_genes = set(patient_df[gene_col].values)
        row = [1 if gene in patient_genes else 0 for gene in top_gene_names]
        matrix.append(row)

    matrix = np.array(matrix).T  # genes = rows, patients = columns

    # -------------------------------------------------
    # 4) Plot heatmap
    # -------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.heatmap(
        matrix,
        yticklabels=top_gene_names,
        cmap=['white', 'darkblue'],
        cbar_kws={'label': 'Present', 'ticks': [0, 1]},
        linewidths=0.5,
        linecolor='gray',
        ax=ax
    )

    # X-axis labels: first 8 characters of the RENAMED IDs
    short_labels = [pid[:8] for pid in sorted_renamed_ids]
    ax.set_xticks(np.arange(len(short_labels)) + 0.5)
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)

    ax.set_xlabel('Patient', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=12, fontweight='bold')
    ax.set_title(
        f'{subtype_name}: Gene Presence/Absence Matrix\n'
        f'Top {top_n} Most Frequent Genes (n={len(sorted_renamed_ids)} patients)',
        fontsize=18, fontweight='bold', pad=20
    )

    # highlight known markers
    km = set(globals().get('known_markers', []))
    for lab in ax.get_yticklabels():
        if lab.get_text() in km:
            lab.set_color('crimson')
            lab.set_fontweight('bold')

    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig(
        f'{output_dir}/{subtype_name}_presence_absence_heatmap.png',
        dpi=300, bbox_inches='tight'
    )
    plt.close()

    print("Saved presence/absence heatmap (CID renamed to FID for sorting/labels)")

def create_metric_distribution_plots(df_summary, subtype_name, output_dir):
    """Create distribution plots for gene metrics"""
    
    # Determine which metrics are available
    has_rank = 'mean_rank' in df_summary.columns
    has_log2fc = 'mean_log2fc' in df_summary.columns
    has_stromal = 'mean_stromal_exclusivity' in df_summary.columns
    
    # Count available plots (frequency is always available)
    n_plots = 1 + sum([has_rank, has_log2fc, has_stromal])
    
    if n_plots == 1:
        fig, axes = plt.subplots(1, 1, figsize=(7, 5))
        axes = [axes]
    elif n_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes = axes.flatten()
    elif n_plots == 3:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
    
    plot_idx = 0
    
    # Frequency distribution (always present)
    axes[plot_idx].hist(df_summary['frequency_pct'], bins=20, 
                    color='steelblue', edgecolor='black', alpha=0.7)
    axes[plot_idx].set_xlabel('Frequency (%)', fontweight='bold')
    axes[plot_idx].set_ylabel('Number of Genes', fontweight='bold')
    axes[plot_idx].set_title('Gene Frequency Distribution', fontweight='bold')
    axes[plot_idx].grid(alpha=0.3)
    
    # Add threshold lines
    for threshold in FREQUENCY_THRESHOLDS.values():
        axes[plot_idx].axvline(threshold * 100, color='red', 
                          linestyle='--', alpha=0.5)
    plot_idx += 1
    
    # Mean rank distribution (if available)
    if has_rank:
        axes[plot_idx].hist(df_summary['mean_rank'], bins=30, 
                        color='coral', edgecolor='black', alpha=0.7)
        axes[plot_idx].set_xlabel('Mean Rank', fontweight='bold')
        axes[plot_idx].set_ylabel('Number of Genes', fontweight='bold')
        axes[plot_idx].set_title('Mean Rank Distribution', fontweight='bold')
        axes[plot_idx].grid(alpha=0.3)
        plot_idx += 1
    
    # Log2FC distribution (if available)
    if has_log2fc:
        axes[plot_idx].hist(df_summary['mean_log2fc'], bins=30, 
                        color='seagreen', edgecolor='black', alpha=0.7)
        axes[plot_idx].set_xlabel('Mean log2FC vs All', fontweight='bold')
        axes[plot_idx].set_ylabel('Number of Genes', fontweight='bold')
        axes[plot_idx].set_title('log2FC Distribution', fontweight='bold')
        axes[plot_idx].grid(alpha=0.3)
        plot_idx += 1
    
    # Stromal exclusivity distribution (if available)
    if has_stromal:
        axes[plot_idx].hist(df_summary['mean_stromal_exclusivity'], bins=30, 
                        color='mediumpurple', edgecolor='black', alpha=0.7)
        axes[plot_idx].set_xlabel('Mean Stromal Exclusivity', fontweight='bold')
        axes[plot_idx].set_ylabel('Number of Genes', fontweight='bold')
        axes[plot_idx].set_title('Stromal Exclusivity Distribution', fontweight='bold')
        axes[plot_idx].grid(alpha=0.3)
        plot_idx += 1
    
    # Hide unused subplots
    for idx in range(plot_idx, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle(f'{subtype_name}: Gene Metric Distributions', 
                 fontsize=16, fontweight='bold', y=1.00)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{subtype_name}_metric_distributions.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved metric distribution plots")

def create_summary_comparison(all_subtypes_data, output_dir):
    """Create comparison plot across all subtypes"""
    
    if len(all_subtypes_data) < 2:
        return
    
    # Prepare data for comparison
    comparison_data = []
    
    for subtype, data in all_subtypes_data.items():
        for category, threshold in FREQUENCY_THRESHOLDS.items():
            n_genes = len(data['categorized'][category])
            comparison_data.append({
                'Subtype': subtype,
                'Category': category,
                'N_Genes': n_genes,
                'Threshold': threshold
            })
    
    df_comparison = pd.DataFrame(comparison_data)
    
    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    categories = list(FREQUENCY_THRESHOLDS.keys())
    x = np.arange(len(all_subtypes_data))
    width = 0.2
    
    for i, category in enumerate(categories):
        category_data = df_comparison[df_comparison['Category'] == category]
        values = [category_data[category_data['Subtype'] == st]['N_Genes'].values[0] 
                 if len(category_data[category_data['Subtype'] == st]) > 0 else 0
                 for st in all_subtypes_data.keys()]
        
        ax.bar(x + i*width, values, width, 
              label=f'{category} (≥{FREQUENCY_THRESHOLDS[category]*100:.0f}%)')
    
    ax.set_xlabel('Cancer Subtype', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('Gene Frequency Comparison Across Subtypes', 
                fontsize=14, fontweight='bold')
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(all_subtypes_data.keys(), rotation=45, ha='right')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/All_Subtypes_Comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved subtype comparison plot")


def create_venn_style_comparison(all_subtypes_data, output_dir):
    """Create Venn-style comparison showing shared genes"""
    
    if len(all_subtypes_data) != 2 and len(all_subtypes_data) != 3:
        print("  Venn comparison only for 2-3 subtypes")
        return
    
    from matplotlib_venn import venn2, venn3
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get universal genes from each subtype
    subtype_names = list(all_subtypes_data.keys())
    gene_sets = [set(all_subtypes_data[st]['categorized']['Frequent']) 
                 for st in subtype_names]
    
    if len(gene_sets) == 2:
        venn2(gene_sets, set_labels=subtype_names, ax=ax)
    elif len(gene_sets) == 3:
        venn3(gene_sets, set_labels=subtype_names, ax=ax)
    
    ax.set_title('Shared Frequent Genes (≥75% frequency) Across Subtypes',
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Subtypes_Venn_Diagram.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved Venn diagram")



################################################################################
# FIG4/FIG5 VISUALS (INLINED): Subtype Exclusivity Heatmap + Primary vs Met Scatter
################################################################################

def create_primary_vs_metastatic_plot(data, output_dir, top_n=25,
                                      primary_label='Primary', metastatic_label='Metastatic',
                                      filename='fig4_primary_vs_metastatic.png'):
    """Figure 4: Primary vs metastatic exclusivity scatter (no category coloring).
    Labels the top N genes by mean exclusivity across the two groups.
    Expects: data['comparison'] with columns:
        - gene
        - exclusivity_mean_primary
        - exclusivity_mean_metastatic
    """
    comparison_df = data.get('comparison')
    if comparison_df is None or comparison_df.empty:
        print("  No comparison data - skipping Fig4")
        return

    df = comparison_df.copy()

    # Coerce numeric columns safely
    df['exclusivity_mean_primary'] = pd.to_numeric(df['exclusivity_mean_primary'], errors='coerce').fillna(0.0)
    df['exclusivity_mean_metastatic'] = pd.to_numeric(df['exclusivity_mean_metastatic'], errors='coerce').fillna(0.0)

    df['mean_exclusivity'] = (df['exclusivity_mean_primary'] + df['exclusivity_mean_metastatic']) / 2.0
    top = df.sort_values('mean_exclusivity', ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(df['exclusivity_mean_primary'],
               df['exclusivity_mean_metastatic'],
               alpha=0.5, s=35, c='gray')

    # Label top genes
    ax.scatter(top['exclusivity_mean_primary'],
               top['exclusivity_mean_metastatic'],
               s=70, c='black', alpha=0.8)

    angles = np.linspace(0, 2*np.pi, len(top), endpoint=False)

    for angle, (_, row) in zip(angles, top.iterrows()):
        dx = 12 * np.cos(angle)
        dy = 12 * np.sin(angle)

        ax.annotate(
            row['gene'],
            xy=(row['exclusivity_mean_primary'], row['exclusivity_mean_metastatic']),
            xytext=(dx, dy),
            textcoords='offset points',
            fontsize=8,
            ha='center',
            va='center',
            arrowprops=dict(
                arrowstyle='-',
                linewidth=0.5,
                color='gray',
                alpha=0.6
            )
        )

    max_val = max(df['exclusivity_mean_primary'].max(),
                  df['exclusivity_mean_metastatic'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)

    ax.set_xlabel(f'{primary_label} Exclusivity', fontweight='bold')
    ax.set_ylabel(f'{metastatic_label} Exclusivity', fontweight='bold')
    ax.set_title(f'{primary_label} vs {metastatic_label}: Stromal Exclusivity', fontweight='bold', fontsize=12)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {filename}")


def create_subtype_comparison_heatmap(data, output_dir, top_n=25,
                                     filename='fig5_subtype_comparison_heatmap.png'):
    """Figure 5: Subtype comparison heatmap (unconstrained color axis).
    Labels top N genes by mean exclusivity across subtypes.
    Expects:
      - data['subtype_comparison'] : dataframe with columns ['gene', <subtype cols>, 'mean_exclusivity']
      - data['available_subtypes'] : list of subtype/group column names to plot
    """
    subtype_comp = data.get('subtype_comparison')
    available_subtypes = data.get('available_subtypes', [])

    if subtype_comp is None or subtype_comp.empty:
        print("  No subtype comparison data - skipping Fig5")
        return
    if not available_subtypes:
        # fall back: infer subtype columns from dataframe
        available_subtypes = [c for c in subtype_comp.columns if c not in ['gene', 'mean_exclusivity']]

     # Top N genes by mean_exclusivity
    if 'mean_exclusivity' in subtype_comp.columns:
        ranked = subtype_comp.sort_values('mean_exclusivity', ascending=False)
    else:
        tmp = subtype_comp.copy()
        cols = [c for c in tmp.columns if c not in ['gene', 'mean_exclusivity']]
        tmp['mean_exclusivity'] = tmp[cols].mean(axis=1)
        ranked = tmp.sort_values('mean_exclusivity', ascending=False)

    top_genes_df = ranked.head(top_n).copy()

    # --- tack on FERMT2 at the end (if present and not already included) ---
    fermt2_mask = ranked['gene'].astype(str).eq('FERMT2')
    if fermt2_mask.any():
        if 'FERMT2' not in top_genes_df['gene'].astype(str).values:
            fermt2_row = ranked.loc[fermt2_mask].head(1)
            top_genes_df = pd.concat([top_genes_df, fermt2_row], ignore_index=True)
    else:
        print("  FERMT2 not found in subtype_comp; not appended.")

    gene_names = top_genes_df['gene'].astype(str).tolist()
    matrix = top_genes_df.set_index('gene')[available_subtypes].copy()

    # Coerce all values numeric
    for c in matrix.columns:
        matrix[c] = pd.to_numeric(matrix[c], errors='coerce').fillna(0.0)

    fig, ax = plt.subplots(figsize=(3, 4))

    # NOTE: unconstrained color axis -> no vmin/vmax/center
    sns.heatmap(matrix,
                cmap='Reds',
                cbar_kws={'label': 'Stromal Exclusivity'},
                ax=ax,
                linewidths=0.5,
                linecolor='white',
                annot=False)

    ax.set_xlabel('Group / Subtype', fontweight='bold')
    ax.set_ylabel('Gene', fontweight='bold')
    ax.set_title('Stromal Exclusivity Across Groups/Subtypes', fontweight='bold', fontsize=12)

    # highlight known markers
    km = set(globals().get('known_markers', []))
    for lab in ax.get_yticklabels():
        lab.set_fontsize(10)
        if lab.get_text() in km:
            lab.set_color('crimson')
            lab.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {filename}")




def create_subtype_exclusivity_swarmplot(data, output_dir,
                                        filename='fig5_subtype_exclusivity_swarmplot.png',
                                        top_n=25,
                                        include_gene='FERMT2'):
    """Swarmplot companion for Fig5.

    Each dot = mean exclusivity of a gene within a subtype/group.
      - x: genes ordered by mean exclusivity (highest -> lowest), with include_gene optionally appended at the end
      - y: mean stromal exclusivity
      - hue: subtype/group

    Expects:
      data['subtype_comparison']: dataframe with columns ['gene', <subtype cols>, 'mean_exclusivity'(optional)]
      data['available_subtypes']: list of subtype columns (optional; inferred if missing)
    """
    subtype_comp = data.get('subtype_comparison')
    available_subtypes = data.get('available_subtypes', [])

    if subtype_comp is None or subtype_comp.empty:
        print("  No subtype comparison data - skipping swarmplot")
        return

    # infer subtype columns if not provided
    if not available_subtypes:
        available_subtypes = [c for c in subtype_comp.columns if c not in ['gene', 'mean_exclusivity']]

    df_all = subtype_comp.copy()

    # Ensure mean_exclusivity exists so we can order genes left->right
    if 'mean_exclusivity' not in df_all.columns:
        df_all['mean_exclusivity'] = (
            df_all[available_subtypes]
            .apply(pd.to_numeric, errors='coerce')
            .fillna(0.0)
            .mean(axis=1)
        )

    ranked = df_all.sort_values('mean_exclusivity', ascending=False).copy()
    top_df = ranked.head(top_n).copy()

    # Append include_gene at the end if present in full table but not in top
    if include_gene is not None:
        mask = ranked['gene'].astype(str).eq(str(include_gene))
        if mask.any() and (str(include_gene) not in top_df['gene'].astype(str).values):
            top_df = pd.concat([top_df, ranked.loc[mask].head(1)], ignore_index=True)

    # Long format for seaborn
    long_df = top_df.melt(
        id_vars=['gene', 'mean_exclusivity'],
        value_vars=available_subtypes,
        var_name='Subtype',
        value_name='Mean_Exclusivity'
    )
    long_df['Mean_Exclusivity'] = pd.to_numeric(long_df['Mean_Exclusivity'], errors='coerce').fillna(0.0)

    # Order genes by mean exclusivity (highest -> lowest), preserving appended include_gene at end
    gene_order = top_df['gene'].astype(str).tolist()

    fig_w = max(10, 0.55 * len(gene_order))
    fig_h = 6
    plt.figure(figsize=(fig_w, fig_h))

    ax = sns.swarmplot(
        data=long_df,
        x='gene',
        y='Mean_Exclusivity',
        hue='Subtype',
        order=gene_order,
        dodge=True,
        size=4,
        alpha=0.9
    )

    ax.set_xlabel("Gene (ranked by mean exclusivity →)", fontweight='bold')
    ax.set_ylabel("Mean Stromal Exclusivity", fontweight='bold')
    ax.set_title("Stromal Exclusivity per Subtype (top genes)", fontweight='bold')

    ax.grid(axis='y', alpha=0.25)
    plt.xticks(rotation=60, ha='right')

    # Legend outside
    ax.legend(title="Subtype", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {filename}")


def build_subtype_comparison_dataframe(all_subtypes_data, metric_col='mean_stromal_exclusivity'):
    """Build a subtype comparison dataframe for Fig5 heatmap/swarmplot.

    Returns a dataframe with columns:
      - gene
      - one column per subtype/group (values = metric_col)
      - mean_exclusivity (mean across subtype/group columns)

    NOTE: This function intentionally returns the FULL gene list (no top-N truncation),
    so downstream plotting code can choose top-N and optionally append specific genes
    like FERMT2 even if they are not in the top set.
    """
    if not all_subtypes_data or len(all_subtypes_data) < 2:
        return None

    subtypes = list(all_subtypes_data.keys())

    # Ensure metric is available in at least one subtype
    if all(metric_col not in all_subtypes_data[st]['summary'].columns for st in subtypes):
        return None

    all_genes = sorted(set().union(*[set(all_subtypes_data[st]['summary']['gene']) for st in subtypes]))
    comp = pd.DataFrame({'gene': all_genes})

    for st in subtypes:
        summ = all_subtypes_data[st]['summary']
        if metric_col in summ.columns:
            gene_to_val = (
                summ[['gene', metric_col]]
                .drop_duplicates('gene')
                .set_index('gene')[metric_col]
            )
            comp[st] = pd.to_numeric(comp['gene'].map(gene_to_val), errors='coerce').fillna(0.0).astype(float)
        else:
            comp[st] = 0.0

    comp['mean_exclusivity'] = comp[subtypes].mean(axis=1)
    comp = comp.sort_values('mean_exclusivity', ascending=False).reset_index(drop=True)
    return comp


def build_primary_vs_metastatic_comparison_dataframe(all_subtypes_data,
                                                     metric_col='mean_stromal_exclusivity'):
    """Build a primary vs metastatic comparison dataframe for Fig4 scatter.
    Only runs when there are exactly 2 groups in all_subtypes_data.
    Returns: (comparison_df, primary_name, metastatic_name)
      comparison_df columns:
        - gene
        - exclusivity_mean_primary
        - exclusivity_mean_metastatic
    """
    if not all_subtypes_data or len(all_subtypes_data) != 2:
        return None

    group_names = list(all_subtypes_data.keys())
    lower = [g.lower() for g in group_names]

    if any('primary' in g for g in lower) and any(('met' in g) or ('metast' in g) for g in lower):
        primary_name = group_names[lower.index(next(g for g in lower if 'primary' in g))]
        metastatic_name = group_names[lower.index(next(g for g in lower if ('met' in g) or ('metast' in g)))]
    else:
        primary_name, metastatic_name = sorted(group_names)

    df_p = all_subtypes_data[primary_name]['summary']
    df_m = all_subtypes_data[metastatic_name]['summary']

    if metric_col not in df_p.columns or metric_col not in df_m.columns:
        return None

    all_genes = sorted(set(df_p['gene']) | set(df_m['gene']))

    p_map = df_p[['gene', metric_col]].drop_duplicates('gene').set_index('gene')[metric_col].to_dict()
    m_map = df_m[['gene', metric_col]].drop_duplicates('gene').set_index('gene')[metric_col].to_dict()

    rows = []
    for g in all_genes:
        rows.append({
            'gene': g,
            'exclusivity_mean_primary': float(p_map.get(g, 0.0) or 0.0),
            'exclusivity_mean_metastatic': float(m_map.get(g, 0.0) or 0.0)
        })

    comp = pd.DataFrame(rows)
    return comp, primary_name, metastatic_name

################################################################################
# MAIN ANALYSIS
################################################################################

################################################################################

def main(base_directory):
    """Main analysis function"""
    
    print("\n" + "="*80)
    print("GENE FREQUENCY ANALYSIS ACROSS CANCER SUBTYPES")
    print("="*80 + "\n")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find subtype folders
    subtype_dirs = find_subtype_folders(base_directory)
    
    if not subtype_dirs:
        print(f"ERROR: No subdirectories found in {base_directory}")
        return
    
    print(f"Found {len(subtype_dirs)} subtype folders:\n")
    
    all_subtypes_data = {}
    
    # Analyze each subtype
    for subtype_dir in subtype_dirs:
        subtype_name = subtype_dir.name
        print(f"{'='*60}")
        print(f"Analyzing: {subtype_name}")
        print(f"{'='*60}")
        
        # Load patient data
        patient_data, patient_ids, gene_col = load_patient_data(subtype_dir)

        if not patient_data or gene_col is None:
            print(f"  No valid CSV files found in {subtype_name}\n")
            continue
        
        n_patients = len(patient_data)
        print(f"  Patients: {n_patients}")
        
        if n_patients < MIN_PATIENTS:
            print(f"  Warning: Only {n_patients} patients (minimum {MIN_PATIENTS} recommended)")
        
        # Analyze gene frequency
        df_summary = analyze_gene_frequency(patient_data, patient_ids, gene_col)
        print(f"  Unique genes: {len(df_summary)}")
        
        # Categorize genes
        print(f"\n  Gene Categories:")
        categorized = categorize_genes(df_summary, n_patients)
        
        # Save summary table
        subtype_output_dir = f"{OUTPUT_DIR}/{subtype_name}"
        os.makedirs(subtype_output_dir, exist_ok=True)
        
        df_summary.to_csv(f'{subtype_output_dir}/gene_frequency_summary.csv', 
                         index=False)
        
        # Save categorized gene lists
        for category, genes in categorized.items():
            with open(f'{subtype_output_dir}/genes_{category}.txt', 'w') as f:
                f.write('\n'.join(genes))
        
        print(f"\n  Creating figures...")
        
        # Create figures
        create_frequency_barplot(df_summary, subtype_name, subtype_output_dir)
        # Additional ranked plots (frequency / mean_rank / composite)
        create_subtype_ranked_plots(df_summary, subtype_name, subtype_output_dir)
        create_rank_vs_frequency_scatter(df_summary, subtype_name, subtype_output_dir)
        create_presence_absence_heatmap(patient_data, patient_ids, 
                                       subtype_name, subtype_output_dir, gene_col)
        create_metric_distribution_plots(df_summary, subtype_name, subtype_output_dir)
        
        # Store for comparison
        all_subtypes_data[subtype_name] = {
            'summary': df_summary,
            'categorized': categorized,
            'n_patients': n_patients
        }
        
        print(f"  Analysis complete for {subtype_name}\n")
    
    # Create comparison figures
    if len(all_subtypes_data) > 1:
        print(f"{'='*60}")
        print("Creating cross-subtype comparisons...")
        print(f"{'='*60}")
        create_summary_comparison(all_subtypes_data, OUTPUT_DIR)

        # Summary plots of top targets across subtypes
        # (1) Shared y-axis genes derived from a global top list
        create_all_subtypes_top_targets_plot(all_subtypes_data, OUTPUT_DIR, top_n=20, rank_metric='composite_score')

        # (2) Per-subtype top lists (each subplot has its own genes)
        create_all_subtypes_top_targets_plot_per_subtype(all_subtypes_data, OUTPUT_DIR, top_n=10, rank_metric='composite_score')
        create_all_subtypes_top_targets_plot_per_subtype(all_subtypes_data, OUTPUT_DIR, top_n=10, rank_metric='frequency_pct')
        create_all_subtypes_top_targets_plot_per_subtype(all_subtypes_data, OUTPUT_DIR, top_n=10, rank_metric='mean_rank')
        
        try:
            create_venn_style_comparison(all_subtypes_data, OUTPUT_DIR)
        except ImportError:
            print("  Install matplotlib-venn for Venn diagrams: pip install matplotlib-venn")
        except Exception as e:
            print(f"  Could not create Venn diagram: {e}")
    
    # ---------------------------------------------------------------------
    # OPTIONAL: Step5-style stromal exclusivity visuals
    #   - Subtype exclusivity heatmap (requires stromal exclusivity metric)
    #   - Primary vs metastatic scatter (ONLY when exactly 2 groups are present)
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # Fig5: Stromal exclusivity heatmap across groups/subtypes (PNG)
    # ---------------------------------------------------------------------
    # Build FULL comparison table (no top-N truncation) so we can always append FERMT2
    comp_df = build_subtype_comparison_dataframe(
        all_subtypes_data,
        metric_col='mean_stromal_exclusivity'
    )
    if comp_df is None:
        print("  Skipping Fig5 subtype plots (missing mean_stromal_exclusivity or <2 groups)")
    else:
        data_for_fig5 = {
            'subtype_comparison': comp_df,
            'available_subtypes': list(all_subtypes_data.keys())
        }
        try:
            create_subtype_comparison_heatmap(
                data_for_fig5,
                OUTPUT_DIR,
                top_n=15,
                filename='fig5_subtype_comparison_heatmap.png'
            )
        except Exception as e:
            print(f"  Could not create Fig5 subtype heatmap: {e}")

        try:
            create_subtype_exclusivity_swarmplot(
                data_for_fig5,
                OUTPUT_DIR,
                filename='fig5_subtype_exclusivity_swarmplot.png',
                top_n=25,
                include_gene='FERMT2'
            )
        except Exception as e:
            print(f"  Could not create Fig5 subtype swarmplot: {e}")

# Fig4: Primary vs metastatic exclusivity scatter (PNG)
    #   Restricted to cases where there are exactly 2 comparison groups
    # ---------------------------------------------------------------------
    if len(all_subtypes_data) == 2:
        comp_res = build_primary_vs_metastatic_comparison_dataframe(
            all_subtypes_data,
            metric_col='mean_stromal_exclusivity'
        )
        if comp_res is None:
            print("  Skipping Fig4 primary vs metastatic plot (missing mean_stromal_exclusivity)")
        else:
            comp_df, primary_name, metastatic_name = comp_res
            data_for_fig4 = {'comparison': comp_df}
            try:
                create_primary_vs_metastatic_plot(
                    data_for_fig4,
                    OUTPUT_DIR,
                    top_n=25,
                    primary_label=primary_name,
                    metastatic_label=metastatic_name,
                    filename='fig4_primary_vs_metastatic.png'
                )
            except Exception as e:
                print(f"  Could not create Fig4 primary vs metastatic plot: {e}")


    # Create final summary report
    print(f"\n{'='*80}")
    print("FINAL SUMMARY")
    print(f"{'='*80}\n")
    
    summary_report = []
    
    for subtype, data in all_subtypes_data.items():
        summary_report.append({
            'Subtype': subtype,
            'N_Patients': data['n_patients'],
            'Total_Genes': len(data['summary']),
            'Universal_100%': len(data['categorized']['Universal']),
            'Frequent_75%': len(data['categorized']['Frequent']),
            'Very_Common_50%': len(data['categorized']['Very_Common']),
            'Common_25%': len(data['categorized']['Common'])
        })
    
    df_final_summary = pd.DataFrame(summary_report)
    df_final_summary.to_csv(f'{OUTPUT_DIR}/FINAL_SUMMARY.csv', index=False)
    
    print(df_final_summary.to_string(index=False))
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved in: {OUTPUT_DIR}/")
    print(f"\nFor each subtype:")
    print(f"  • gene_frequency_summary.csv - Complete gene statistics")
    print(f"  • genes_[Category].txt - Gene lists by frequency category")
    print(f"  • [Subtype]_frequency_barplot.png - Top genes by frequency")
    print(f"  • [Subtype]_rank_vs_frequency.png - Rank vs frequency scatter")
    print(f"  • [Subtype]_presence_absence_heatmap.png - Patient-gene matrix")
    print(f"  • [Subtype]_metric_distributions.png - Statistical distributions")
    print(f"\nComparison files:")
    print(f"  • All_Subtypes_Comparison.png - Bar chart comparison")
    print(f"  • FINAL_SUMMARY.csv - Summary table\n")

################################################################################
# RUN ANALYSIS
################################################################################

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        base_directory = sys.argv[1]
    else:
        # Prompt for directory
        base_directory = input("Enter the path to the base directory containing subtype folders: ").strip()
    
    if not os.path.exists(base_directory):
        print(f"ERROR: Directory not found: {base_directory}")
        sys.exit(1)
    
    main(base_directory)