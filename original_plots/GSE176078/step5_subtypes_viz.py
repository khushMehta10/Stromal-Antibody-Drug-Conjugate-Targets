"""
STEP 5: Create Visualizations for Subtype Analysis
===================================================
Creates publication-quality figures for breast cancer subtype comparison

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
        'ER': 'ER_stromal_targets.csv',
        'HER2': 'HER2_stromal_targets.csv',
        'TNBC': 'TNBC_stromal_targets.csv',
        'comparison': 'subtype_comparison.csv'
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
    ax.axhline(20, color='blue', linestyle='--', alpha=0.5, label='Exclusivity threshold')
    ax.axvline(5.0, color='green', linestyle='--', alpha=0.5, label='Log2FC threshold')
    
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
    Figure 3: Heatmap of top targets across patients (colored by subtype)
    """
    print("\nCreating Figure 3: Patient heatmap...")
    
    per_patient_df = data['per_patient']
    patient_stats = data['patient_stats']
    
    if per_patient_df is None:
        return
    
    # Get top 20 genes from global
    if data['global'] is not None:
        top_genes = data['global'].head(20)['gene'].tolist()
    else:
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
    
    # Sort columns by subtype if available
    if patient_stats is not None and 'sample_type' in patient_stats.columns:
        # Get patient subtypes
        patient_subtypes = patient_stats.set_index('patient_id')['sample_type']
        
        # Sort patients by subtype
        patient_order = []
        for subtype in ['ER+', 'HER2+', 'TNBC']:
            subtype_patients = patient_subtypes[patient_subtypes == subtype].index.tolist()
            patient_order.extend([p for p in subtype_patients if p in matrix.columns])
        
        # Add any remaining patients
        for p in matrix.columns:
            if p not in patient_order:
                patient_order.append(p)
        
        matrix = matrix[patient_order]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(16, 8))
    
    sns.heatmap(matrix, cmap='YlOrRd', center=0, 
               cbar_kws={'label': 'Stromal Exclusivity'},
               ax=ax, linewidths=0.5, linecolor='white')
    
    # Add subtype color bar if available
    if patient_stats is not None and 'sample_type' in patient_stats.columns:
        subtype_colors = {'ER+': 'lightblue', 'HER2+': 'pink', 'TNBC': 'lightgreen'}
        
        # Create colored boxes above heatmap
        for i, patient in enumerate(matrix.columns):
            if patient in patient_subtypes.index:
                subtype = patient_subtypes[patient]
                color = subtype_colors.get(subtype, 'gray')
                ax.add_patch(plt.Rectangle((i, -0.5), 1, 0.5, 
                                          facecolor=color, edgecolor='black', 
                                          clip_on=False))
    
    ax.set_xlabel('Patient ID (colored by subtype)', fontweight='bold')
    ax.set_ylabel('Gene', fontweight='bold')
    ax.set_title('Stromal Target Exclusivity Across Patients', fontweight='bold', fontsize=12)
    
    # Add legend for subtypes
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', label='ER+'),
        Patch(facecolor='pink', label='HER2+'),
        Patch(facecolor='lightgreen', label='TNBC')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig3_patient_heatmap.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig3_patient_heatmap.pdf")


def create_subtype_comparison_plot(data, output_dir):
    """
    Figure 4: Subtype comparison (3-way plot or heatmap)
    """
    print("\nCreating Figure 4: Subtype comparison...")
    
    comparison_df = data['comparison']
    if comparison_df is None:
        return
    
    # Get top 15 genes by average exclusivity
    top_genes = comparison_df.head(15).copy()
    
    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    x = np.arange(len(top_genes))
    width = 0.25
    
    bars1 = ax.bar(x - width, top_genes['exclusivity_mean_ER'], width, 
                   label='ER+', color='lightblue', edgecolor='black')
    bars2 = ax.bar(x, top_genes['exclusivity_mean_HER2'], width, 
                   label='HER2+', color='pink', edgecolor='black')
    bars3 = ax.bar(x + width, top_genes['exclusivity_mean_TNBC'], width, 
                   label='TNBC', color='lightgreen', edgecolor='black')
    
    ax.set_xlabel('Gene', fontweight='bold')
    ax.set_ylabel('Stromal Exclusivity Score', fontweight='bold')
    ax.set_title('Top 15 Stromal Targets Across Breast Cancer Subtypes', 
                 fontweight='bold', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(top_genes['gene'], rotation=45, ha='right')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig4_subtype_comparison.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig4_subtype_comparison.pdf")


def create_subtype_scatter_plot(data, output_dir):
    """
    Figure 5: 3D scatter or pairwise scatter plots
    """
    print("\nCreating Figure 5: Subtype scatter plots...")
    
    comparison_df = data['comparison']
    if comparison_df is None:
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: ER+ vs HER2+
    ax = axes[0]
    ax.scatter(comparison_df['exclusivity_mean_ER'],
              comparison_df['exclusivity_mean_HER2'],
              alpha=0.5, s=50, c='gray')
    
    # Highlight universal
    universal = comparison_df[comparison_df['target_type'] == 'universal']
    if len(universal) > 0:
        ax.scatter(universal['exclusivity_mean_ER'],
                  universal['exclusivity_mean_HER2'],
                  s=100, c='purple', edgecolors='black', linewidth=1,
                  label='Universal', zorder=5)
        
        # Label top 5
        for _, row in universal.head(5).iterrows():
            ax.annotate(row['gene'],
                       (row['exclusivity_mean_ER'], row['exclusivity_mean_HER2']),
                       fontsize=8, ha='left', va='bottom')
    
    max_val = max(comparison_df['exclusivity_mean_ER'].max(),
                  comparison_df['exclusivity_mean_HER2'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)
    ax.set_xlabel('ER+ Exclusivity', fontweight='bold')
    ax.set_ylabel('HER2+ Exclusivity', fontweight='bold')
    ax.set_title('ER+ vs HER2+', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 2: ER+ vs TNBC
    ax = axes[1]
    ax.scatter(comparison_df['exclusivity_mean_ER'],
              comparison_df['exclusivity_mean_TNBC'],
              alpha=0.5, s=50, c='gray')
    
    if len(universal) > 0:
        ax.scatter(universal['exclusivity_mean_ER'],
                  universal['exclusivity_mean_TNBC'],
                  s=100, c='purple', edgecolors='black', linewidth=1,
                  label='Universal', zorder=5)
    
    max_val = max(comparison_df['exclusivity_mean_ER'].max(),
                  comparison_df['exclusivity_mean_TNBC'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)
    ax.set_xlabel('ER+ Exclusivity', fontweight='bold')
    ax.set_ylabel('TNBC Exclusivity', fontweight='bold')
    ax.set_title('ER+ vs TNBC', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 3: HER2+ vs TNBC
    ax = axes[2]
    ax.scatter(comparison_df['exclusivity_mean_HER2'],
              comparison_df['exclusivity_mean_TNBC'],
              alpha=0.5, s=50, c='gray')
    
    if len(universal) > 0:
        ax.scatter(universal['exclusivity_mean_HER2'],
                  universal['exclusivity_mean_TNBC'],
                  s=100, c='purple', edgecolors='black', linewidth=1,
                  label='Universal', zorder=5)
    
    max_val = max(comparison_df['exclusivity_mean_HER2'].max(),
                  comparison_df['exclusivity_mean_TNBC'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)
    ax.set_xlabel('HER2+ Exclusivity', fontweight='bold')
    ax.set_ylabel('TNBC Exclusivity', fontweight='bold')
    ax.set_title('HER2+ vs TNBC', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig5_subtype_scatter.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig5_subtype_scatter.pdf")


def create_summary_statistics_plot(data, output_dir):
    """
    Figure 6: Summary statistics by subtype
    """
    print("\nCreating Figure 6: Summary statistics...")
    
    patient_stats = data['patient_stats']
    if patient_stats is None:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Stromal cells by subtype
    ax = axes[0, 0]
    if 'sample_type' in patient_stats.columns:
        subtype_stromal = patient_stats.groupby('sample_type')['n_stromal_cells'].mean()
        colors = {'ER+': 'lightblue', 'HER2+': 'pink', 'TNBC': 'lightgreen'}
        bar_colors = [colors.get(st, 'gray') for st in subtype_stromal.index]
        
        ax.bar(subtype_stromal.index, subtype_stromal.values, 
               color=bar_colors, edgecolor='black')
        ax.set_ylabel('Mean Stromal Cells per Patient', fontweight='bold')
        ax.set_title('Stromal Cell Abundance by Subtype', fontweight='bold')
    
    # 2. Tumor cells by subtype
    ax = axes[0, 1]
    if 'sample_type' in patient_stats.columns:
        subtype_tumor = patient_stats.groupby('sample_type')['n_tumor_cells'].mean()
        bar_colors = [colors.get(st, 'gray') for st in subtype_tumor.index]
        
        ax.bar(subtype_tumor.index, subtype_tumor.values, 
               color=bar_colors, edgecolor='black')
        ax.set_ylabel('Mean Tumor Cells per Patient', fontweight='bold')
        ax.set_title('Tumor Cell Abundance by Subtype', fontweight='bold')
    
    # 3. Stromal percentage by subtype
    ax = axes[1, 0]
    if 'sample_type' in patient_stats.columns:
        subtype_pct = patient_stats.groupby('sample_type')['pct_stromal'].mean()
        bar_colors = [colors.get(st, 'gray') for st in subtype_pct.index]
        
        ax.bar(subtype_pct.index, subtype_pct.values, 
               color=bar_colors, edgecolor='black')
        ax.set_ylabel('Mean Stromal Percentage (%)', fontweight='bold')
        ax.set_title('Stromal Content by Subtype', fontweight='bold')
    
    # 4. Patient counts
    ax = axes[1, 1]
    if 'sample_type' in patient_stats.columns:
        sample_counts = patient_stats['sample_type'].value_counts()
        bar_colors = [colors.get(st, 'gray') for st in sample_counts.index]
        
        ax.bar(sample_counts.index, sample_counts.values, 
               color=bar_colors, edgecolor='black')
        ax.set_ylabel('Number of Patients', fontweight='bold')
        ax.set_title('Patient Distribution by Subtype', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig6_summary_statistics.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig6_summary_statistics.pdf")


def create_summary_tables(data, output_dir):
    """
    Create summary tables for the paper
    """
    print("\nCreating summary tables...")
    
    # Table 1: Top 30 global targets
    if data['global'] is not None:
        top_targets = data['global'].head(30)[['gene', 'stromal_exclusivity', 'log2fc_vs_all',
                                               'stromal_detection', 'tumor_penalty']].copy()
        
        top_targets['stromal_exclusivity'] = top_targets['stromal_exclusivity'].round(2)
        top_targets['log2fc_vs_all'] = top_targets['log2fc_vs_all'].round(2)
        top_targets['stromal_detection'] = (top_targets['stromal_detection'] * 100).round(1)
        top_targets['tumor_penalty'] = top_targets['tumor_penalty'].round(2)
        
        top_targets.columns = ['Gene', 'Exclusivity', 'Log2FC', 'Detection (%)', 'Tumor Penalty']
        
        top_targets.to_csv(output_dir / 'table1_top_targets.csv', index=False)
        print("  Saved: table1_top_targets.csv")
    
    # Table 2: Subtype-specific targets
    if data['comparison'] is not None:
        comparison = data['comparison']
        
        # Universal targets
        universal = comparison[comparison['target_type'] == 'universal'][
            ['gene', 'exclusivity_mean_ER', 'exclusivity_mean_HER2', 
             'exclusivity_mean_TNBC', 'exclusivity_cv']
        ].head(20).copy()
        
        if len(universal) > 0:
            universal = universal.round(2)
            universal.columns = ['Gene', 'ER+ Exclusivity', 'HER2+ Exclusivity', 
                               'TNBC Exclusivity', 'CV']
            universal.to_csv(output_dir / 'table2_universal_targets.csv', index=False)
            print("  Saved: table2_universal_targets.csv")


def main():
    print("\n" + "="*80)
    print("STEP 5: CREATE VISUALIZATIONS (SUBTYPE ANALYSIS)")
    print("="*80)
    
    # Load results
    data = load_all_results()
    
    # Create output directory
    output_dir = Path('./results/figures')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Create figures
    create_top_targets_barplot(data, output_dir)
    create_volcano_plot(data, output_dir)
    create_patient_heatmap(data, output_dir)
    create_subtype_comparison_plot(data, output_dir)
    create_subtype_scatter_plot(data, output_dir)
    create_summary_statistics_plot(data, output_dir)
    create_summary_tables(data, output_dir)
    
    print("\n" + "="*80)
    print("STEP 5 COMPLETE")
    print("="*80)
    print(f"\nAll figures saved in: {output_dir.absolute()}")
    print("\nSubtype analysis pipeline complete.")
    print("Review the results in ./results/")


if __name__ == "__main__":
    main()
