"""
STEP 5: Create Visualizations
==============================
Creates publication-quality figures for conference presentation

Input: All results from previous steps
Output: ./results/figures/*.pdf

Runtime: ~5 minutes
Memory: ~500MB
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

def load_all_results():
    """Load all result files"""
    results_dir = Path('./results')
    
    data = {}
    
    # Load files
    files = {
        'global': 'global_stromal_targets.csv',
        'per_patient': 'per_patient_stromal_targets.csv',
        'patient_stats': 'patient_statistics.csv',
        'universal': 'universal_stromal_targets.csv',
        'primary': 'primary_stromal_targets.csv',
        'metastatic': 'metastatic_stromal_targets.csv',
        'comparison': 'primary_vs_metastatic_comparison.csv'
    }
    
    for key, filename in files.items():
        filepath = results_dir / filename
        if filepath.exists():
            data[key] = pd.read_csv(filepath)
            print(f"Loaded {filename}")
        else:
            data[key] = None
            print(f"{filename} not found")
    
    return data


def create_top_targets_barplot(data, output_dir):
    """
    Figure 1: Top 20 global stromal targets
    """
    print("\nCreating Figure 1: Top targets barplot...")
    
    global_df = data['global']
    if global_df is None:
        return
    
    top_20 = global_df.head(20).copy()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create barplot
    y_pos = np.arange(len(top_20))
    bars = ax.barh(y_pos, top_20['stromal_exclusivity'], color='steelblue')
    
    # Color known markers differently
    known_markers = ['FAP', 'PDPN', 'COL1A1', 'COL1A2', 'COL11A1', 'POSTN', 'THY1',
                     'PDGFRA', 'PDGFRB', 'VIM', 'DCN', 'LUM', 'ACTA2', 'FN1',
                     'S100A4', 'SPARC', 'TAGLN', 'MMP2', 'MMP11', 'VCAN', 'CDH11']
    for i, gene in enumerate(top_20['gene']):
        if gene in known_markers:
            bars[i].set_color('coral')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_20['gene'])
    ax.invert_yaxis()
    ax.set_xlabel('Stromal Exclusivity Score', fontweight='bold')
    ax.set_title('Top 20 Stromal-Specific ADC Targets', fontweight='bold', fontsize=12)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='coral', label='Known stromal markers'),
        Patch(facecolor='steelblue', label='Novel targets')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig1_top_targets.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig1_top_targets.pdf")


def create_volcano_plot(data, output_dir):
    """
    Figure 2: Volcano plot of log2FC vs exclusivity
    """
    print("\nCreating Figure 2: Volcano plot...")
    
    global_df = data['global']
    if global_df is None:
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot all points
    ax.scatter(global_df['log2fc_vs_all'], 
              global_df['stromal_exclusivity'],
              alpha=0.3, s=20, c='gray')
    
    # Highlight top 20
    top_20 = global_df.head(20)
    ax.scatter(top_20['log2fc_vs_all'],
              top_20['stromal_exclusivity'],
              s=80, c='red', edgecolors='black', linewidth=1,
              label='Top 20 targets', zorder=5)
    
    # Label top 10
    for _, row in top_20.head(10).iterrows():
        ax.annotate(row['gene'], 
                   (row['log2fc_vs_all'], row['stromal_exclusivity']),
                   fontsize=8, ha='left', va='bottom')
    
    # Add threshold lines
    ax.axhline(0.5, color='blue', linestyle='--', alpha=0.5, label='Exclusivity threshold')
    ax.axvline(1.0, color='green', linestyle='--', alpha=0.5, label='Log2FC threshold')
    
    ax.set_xlabel('Log2 Fold Change (Stromal vs Non-Stromal)', fontweight='bold')
    ax.set_ylabel('Stromal Exclusivity Score', fontweight='bold')
    ax.set_title('Stromal Specificity Landscape', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig2_volcano.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig2_volcano.pdf")


def create_patient_heatmap(data, output_dir):
    """
    Figure 3: Heatmap of top targets across patients
    """
    print("\nCreating Figure 3: Patient heatmap...")
    
    per_patient_df = data['per_patient']
    if per_patient_df is None:
        return
    
    # Get top 20 genes from global
    if data['global'] is not None:
        top_genes = data['global'].head(20)['gene'].tolist()
    else:
        # Use most common genes
        gene_counts = per_patient_df['gene'].value_counts()
        top_genes = gene_counts.head(20).index.tolist()
    
    # Filter to top genes
    df_subset = per_patient_df[per_patient_df['gene'].isin(top_genes)]
    
    # Pivot to matrix
    matrix = df_subset.pivot_table(
        index='gene', 
        columns='patient_id', 
        values='stromal_exclusivity',
        fill_value=0
    )
    
    # Reorder genes by mean exclusivity
    gene_order = matrix.mean(axis=1).sort_values(ascending=False).index
    matrix = matrix.loc[gene_order]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(14, 8))
    
    sns.heatmap(matrix, cmap='YlOrRd', center=0, 
               cbar_kws={'label': 'Stromal Exclusivity'},
               ax=ax, linewidths=0.5, linecolor='white')
    
    ax.set_xlabel('Patient ID', fontweight='bold')
    ax.set_ylabel('Gene', fontweight='bold')
    ax.set_title('Stromal Target Exclusivity Across Patients', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig3_patient_heatmap.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig3_patient_heatmap.pdf")


def create_primary_vs_metastatic_plot(data, output_dir):
    """
    Figure 4: Primary vs metastatic comparison
    """
    print("\nCreating Figure 4: Primary vs metastatic...")
    
    comparison_df = data['comparison']
    if comparison_df is None:
        return
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Scatter plot
    ax.scatter(comparison_df['exclusivity_mean_primary'],
              comparison_df['exclusivity_mean_metastatic'],
              alpha=0.5, s=50, c='gray')
    
    # Highlight universal targets
    universal = comparison_df[comparison_df['target_type'] == 'universal']
    ax.scatter(universal['exclusivity_mean_primary'],
              universal['exclusivity_mean_metastatic'],
              s=100, c='green', edgecolors='black', linewidth=1,
              label='Universal', zorder=5)
    
    # Highlight stage-specific
    primary_spec = comparison_df[comparison_df['target_type'] == 'primary_specific']
    if len(primary_spec) > 0:
        ax.scatter(primary_spec['exclusivity_mean_primary'],
                  primary_spec['exclusivity_mean_metastatic'],
                  s=100, c='blue', edgecolors='black', linewidth=1,
                  label='Primary-specific', zorder=5)
    
    met_spec = comparison_df[comparison_df['target_type'] == 'metastatic_specific']
    if len(met_spec) > 0:
        ax.scatter(met_spec['exclusivity_mean_primary'],
                  met_spec['exclusivity_mean_metastatic'],
                  s=100, c='red', edgecolors='black', linewidth=1,
                  label='Metastatic-specific', zorder=5)
    
    # Label top universal targets
    for _, row in universal.head(5).iterrows():
        ax.annotate(row['gene'],
                   (row['exclusivity_mean_primary'], row['exclusivity_mean_metastatic']),
                   fontsize=8, ha='left', va='bottom')
    
    # Diagonal line
    max_val = max(comparison_df['exclusivity_mean_primary'].max(),
                  comparison_df['exclusivity_mean_metastatic'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Equal exclusivity')
    
    ax.set_xlabel('Primary Tumor Exclusivity', fontweight='bold')
    ax.set_ylabel('Metastatic Tumor Exclusivity', fontweight='bold')
    ax.set_title('Primary vs Metastatic Stromal Targets', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig4_primary_vs_metastatic.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig4_primary_vs_metastatic.pdf")


def create_summary_statistics_plot(data, output_dir):
    """
    Figure 5: Summary statistics
    """
    print("\nCreating Figure 5: Summary statistics...")
    
    patient_stats = data['patient_stats']
    if patient_stats is None:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Cells per patient
    ax = axes[0, 0]
    ax.bar(['Stromal', 'Tumor', 'Other'],
          [patient_stats['n_stromal_cells'].mean(),
           patient_stats['n_tumor_cells'].mean(),
           patient_stats['n_total_cells'].mean() - 
           patient_stats['n_stromal_cells'].mean() - 
           patient_stats['n_tumor_cells'].mean()],
          color=['steelblue', 'coral', 'lightgray'])
    ax.set_ylabel('Mean Cells per Patient', fontweight='bold')
    ax.set_title('Cell Type Distribution', fontweight='bold')
    
    # 2. Percentage distribution
    ax = axes[0, 1]
    ax.bar(['Stromal', 'Tumor'],
          [patient_stats['pct_stromal'].mean(),
           patient_stats['pct_tumor'].mean()],
          color=['steelblue', 'coral'])
    ax.set_ylabel('Percentage (%)', fontweight='bold')
    ax.set_title('Cell Type Percentages', fontweight='bold')
    
    # 3. Patients by sample type
    ax = axes[1, 0]
    if 'sample_type' in patient_stats.columns:
        sample_counts = patient_stats['sample_type'].value_counts()
        ax.bar(sample_counts.index, sample_counts.values, color='mediumpurple')
        ax.set_ylabel('Number of Patients', fontweight='bold')
        ax.set_title('Patients by Sample Type', fontweight='bold')
    
    # 4. Detection rate distribution
    ax = axes[1, 1]
    if data['global'] is not None:
        ax.hist(data['global']['stromal_detection'] * 100, 
               bins=30, color='steelblue', edgecolor='black')
        ax.set_xlabel('Stromal Detection Rate (%)', fontweight='bold')
        ax.set_ylabel('Number of Genes', fontweight='bold')
        ax.set_title('Detection Rate Distribution', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig5_summary_statistics.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig5_summary_statistics.pdf")


def create_summary_table(data, output_dir):
    """
    Create a summary table for the paper
    """
    print("\nCreating summary table...")
    
    if data['global'] is None:
        return
    
    # Get top 30 targets
    top_targets = data['global'].head(30)[['gene', 'stromal_exclusivity', 'log2fc_vs_all',
                                           'stromal_detection', 'tumor_penalty']].copy()
    
    # Round values
    top_targets['stromal_exclusivity'] = top_targets['stromal_exclusivity'].round(3)
    top_targets['log2fc_vs_all'] = top_targets['log2fc_vs_all'].round(2)
    top_targets['stromal_detection'] = (top_targets['stromal_detection'] * 100).round(1)
    top_targets['tumor_penalty'] = top_targets['tumor_penalty'].round(2)
    
    # Rename columns
    top_targets.columns = ['Gene', 'Exclusivity', 'Log2FC', 'Detection (%)', 'Tumor Penalty']
    
    # Save as CSV
    top_targets.to_csv(output_dir / 'table1_top_targets.csv', index=False)
    
    print("  Saved: table1_top_targets.csv")


def main():
    print("\n" + "="*80)
    print("STEP 5: CREATE VISUALIZATIONS")
    print("="*80)
    
    # Load results
    data = load_all_results()
    
    # Create output directory
    output_dir = Path('./results/figures')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Create figures
    create_top_targets_barplot(data, output_dir)
    create_volcano_plot(data, output_dir)
    # create_patient_heatmap(data, output_dir)
    create_primary_vs_metastatic_plot(data, output_dir)
    create_summary_statistics_plot(data, output_dir)
    create_summary_table(data, output_dir)
    
    print("\n" + "="*80)
    print("STEP 5 COMPLETE")
    print("="*80)
    print(f"\nAll figures saved in: {output_dir.absolute()}")
    print("\nAnalysis pipeline complete.")
    print("Review the results in ./results/")


if __name__ == "__main__":
    main()
