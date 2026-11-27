"""
STEP 5: Create Visualizations - SUBTYPE-AWARE
==============================================
Creates publication-quality figures including subtype analysis

Input: All results from previous steps
Output: ./results/figures/*.pdf

Runtime: ~5-10 minutes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for publication

def load_all_results():
    """Load all result files including subtype-specific ones"""
    results_dir = Path('./results')
    
    data = {}
    
    # Core files
    files = {
        'global': 'global_stromal_targets.csv',
        'per_patient': 'per_patient_stromal_targets.csv',
        'patient_stats': 'patient_statistics.csv',
        'universal': 'universal_stromal_targets.csv',
        'primary': 'primary_stromal_targets.csv',
        'metastatic': 'metastatic_stromal_targets.csv',
        'comparison': 'primary_vs_metastatic_comparison.csv'
    }
    
    # Subtype-specific files
    subtype_files = {
        'subtype_comparison': 'subtype_comparison.csv',
        'per_patient_by_subtype': 'per_patient_by_subtype.csv',
        'primary_vs_met_by_subtype': 'primary_vs_metastatic_by_subtype.csv',
        'ER+': 'ER+_stromal_targets.csv',
        'TNBC': 'TNBC_stromal_targets.csv',
        'HER2+': 'HER2+_stromal_targets.csv',
        'Normal': 'Normal_stromal_targets.csv'
    }
    
    all_files = {**files, **subtype_files}
    
    for key, filename in all_files.items():
        filepath = results_dir / filename
        if filepath.exists():
            data[key] = pd.read_csv(filepath)
            print(f"Loaded {filename}")
        else:
            data[key] = None
            if key in files:  # Core files
                print(f"{filename} not found")
    
    # Detect available subtypes
    available_subtypes = [k for k in ['ER+', 'TNBC', 'HER2+', 'Normal'] if data.get(k) is not None]
    data['available_subtypes'] = available_subtypes
    
    if available_subtypes:
        print(f"\nSubtype data found for: {', '.join(available_subtypes)}")
    else:
        print(f"\nNo subtype-specific data found")
    
    return data


def create_top_targets_barplot(data, output_dir):
    """Figure 1: Top 20 global targets"""
    print("\nCreating Figure 1: Top global targets...")
    
    global_df = data['global']
    if global_df is None:
        return
    
    top_20 = global_df.head(20).copy()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(top_20))
    bars = ax.barh(y_pos, top_20['stromal_exclusivity'], color='steelblue')
    
    # Known stromal/CAF markers - comprehensive list
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
    ax.set_title('Top 20 Global Stromal-Specific ADC Targets', fontweight='bold', fontsize=12)
    
    legend_elements = [
        Patch(facecolor='coral', label='Known stromal markers'),
        Patch(facecolor='steelblue', label='Novel targets')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig1_top_global_targets.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig1_top_global_targets.pdf")


def create_volcano_plot(data, output_dir):
    """Figure 2: Volcano plot"""
    print("\nCreating Figure 2: Volcano plot...")
    
    global_df = data['global']
    if global_df is None:
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.scatter(global_df['log2fc_vs_all'], 
              global_df['stromal_exclusivity'],
              alpha=0.3, s=20, c='gray', label='All genes')
    
    top_20 = global_df.head(20)
    ax.scatter(top_20['log2fc_vs_all'],
              top_20['stromal_exclusivity'],
              s=80, c='red', edgecolors='black', linewidth=1,
              label='Top 20 targets', zorder=5)
    
    for _, row in top_20.head(10).iterrows():
        ax.annotate(row['gene'], 
                   (row['log2fc_vs_all'], row['stromal_exclusivity']),
                   fontsize=8, ha='left', va='bottom')
    
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
    """Figure 3: Patient heatmap with subtype annotation"""
    print("\nCreating Figure 3: Patient heatmap...")
    
    per_patient_df = data['per_patient']
    patient_stats = data['patient_stats']
    
    if per_patient_df is None:
        return
    
    # Get top genes
    if data['global'] is not None:
        top_genes = data['global'].head(20)['gene'].tolist()
    else:
        gene_counts = per_patient_df['gene'].value_counts()
        top_genes = gene_counts.head(20).index.tolist()
    
    df_subset = per_patient_df[per_patient_df['gene'].isin(top_genes)]
    
    matrix = df_subset.pivot_table(
        index='gene', 
        columns='patient_id', 
        values='stromal_exclusivity',
        fill_value=0
    )
    
    gene_order = matrix.mean(axis=1).sort_values(ascending=False).index
    matrix = matrix.loc[gene_order]
    
    # Create figure with subtype annotation if available
    if patient_stats is not None and 'cancer_subtype' in patient_stats.columns:
        fig = plt.figure(figsize=(16, 8))
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 15], hspace=0.02)
        
        ax_subtype = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1])
        
        # Subtype annotation with numeric mapping
        patient_subtypes = patient_stats.set_index('patient_id')['cancer_subtype']
        
        # Map subtypes to numbers for imshow
        subtype_to_num = {'ER+': 0, 'TNBC': 1, 'HER2+': 2, 'Normal': 3, 'Unknown': 4}
        subtype_to_color = {0: 'pink', 1: 'purple', 2: 'orange', 3: 'lightgreen', 4: 'gray'}
        print(matrix)
        # Match matrix column order and convert to numbers
        ordered_subtypes = [patient_subtypes.get(p, 'Unknown') for p in matrix.columns]
        numeric_subtypes = [subtype_to_num.get(s, 4) for s in ordered_subtypes]
        
        # Create custom colormap
        from matplotlib.colors import ListedColormap
        colors = [subtype_to_color[i] for i in range(5)]
        cmap_subtypes = ListedColormap(colors)
        
        # Plot subtype bar
        ax_subtype.imshow([numeric_subtypes], aspect='auto', cmap=cmap_subtypes, 
                         vmin=0, vmax=4, interpolation='nearest')
        ax_subtype.set_xticks([])
        ax_subtype.set_yticks([])
        ax_subtype.set_ylabel('Subtype', fontsize=8)
        
        # Legend for subtypes - only show subtypes that exist
        unique_subtypes = set(ordered_subtypes)
        legend_elements = [Patch(facecolor=subtype_to_color[subtype_to_num[s]], label=s) 
                          for s in ['ER+', 'TNBC', 'HER2+', 'Normal'] 
                          if s in unique_subtypes]
        ax_subtype.legend(handles=legend_elements, loc='center left', 
                         bbox_to_anchor=(1, 0.5), frameon=False, fontsize=8)
        
        # Heatmap
        sns.heatmap(matrix, cmap='YlOrRd', center=0, 
                   cbar_kws={'label': 'Stromal Exclusivity'},
                   ax=ax_heatmap, linewidths=0.5, linecolor='white',
                   xticklabels=True, yticklabels=True)
        
        ax_heatmap.set_xlabel('Patient ID', fontweight='bold')
        ax_heatmap.set_ylabel('Gene', fontweight='bold')
        
        fig.suptitle('Stromal Target Exclusivity Across Patients (by Subtype)', 
                    fontweight='bold', fontsize=12, y=0.98)
    else:
        # Simple heatmap without subtype
        fig, ax = plt.subplots(figsize=(14, 8))
        
        sns.heatmap(matrix, cmap='YlOrRd', center=0, 
                   cbar_kws={'label': 'Stromal Exclusivity'},
                   ax=ax, linewidths=0.5, linecolor='white')
        
        ax.set_xlabel('Patient ID', fontweight='bold')
        ax.set_ylabel('Gene', fontweight='bold')
        ax.set_title('Stromal Target Exclusivity Across Patients', 
                    fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig3_patient_heatmap.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig3_patient_heatmap.pdf")


def create_primary_subtype_plot(data, output_dir):
    """Alternative Figure 4: Compare primary tumor subtypes"""
    available_subtypes = data['available_subtypes']
    if not available_subtypes or len(available_subtypes) < 2:
        print("  Need at least 2 subtypes for comparison")
        return
    
    # Get data for each subtype
    subtype_data = {}
    for subtype in available_subtypes:
        if data[subtype] is not None:
            top_20 = data[subtype].head(20)['gene'].tolist()
            subtype_data[subtype] = {gene: data[subtype][data[subtype]['gene'] == gene]['stromal_exclusivity'].values[0] 
                                    for gene in top_20}
    
    if len(subtype_data) < 2:
        return
    
    # Create scatter plot comparing two subtypes
    subtypes = list(subtype_data.keys())
    subtype1, subtype2 = subtypes[0], subtypes[1]
    
    # Get common genes
    all_genes = set(subtype_data[subtype1].keys()) | set(subtype_data[subtype2].keys())
    
    x_vals = [subtype_data[subtype1].get(gene, 0) for gene in all_genes]
    y_vals = [subtype_data[subtype2].get(gene, 0) for gene in all_genes]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    ax.scatter(x_vals, y_vals, alpha=0.6, s=80, c='steelblue', edgecolors='black', linewidth=0.5)
    
    # Annotate top genes
    for gene in list(all_genes)[:10]:
        x = subtype_data[subtype1].get(gene, 0)
        y = subtype_data[subtype2].get(gene, 0)
        if x > 0.5 or y > 0.5:
            ax.annotate(gene, (x, y), fontsize=8, ha='left', va='bottom')
    
    max_val = max(max(x_vals), max(y_vals))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Equal exclusivity')
    
    ax.set_xlabel(f'{subtype1} Stromal Exclusivity', fontweight='bold')
    ax.set_ylabel(f'{subtype2} Stromal Exclusivity', fontweight='bold')
    ax.set_title(f'Primary Tumor Stromal Targets: {subtype1} vs {subtype2}', 
                fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig4_primary_subtype_comparison.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig4_primary_subtype_comparison.pdf")


def create_primary_vs_metastatic_plot(data, output_dir):
    """Figure 4: Primary vs metastatic (or primary subtypes if no metastatic)"""
    print("\nCreating Figure 4: Primary analysis...")
    
    comparison_df = data['comparison']
    
    # Check if we have metastatic comparison
    if comparison_df is None:
        print("  No metastatic data - creating primary subtype comparison instead")
        create_primary_subtype_plot(data, output_dir)
        return
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    ax.scatter(comparison_df['exclusivity_mean_primary'],
              comparison_df['exclusivity_mean_metastatic'],
              alpha=0.5, s=50, c='gray', label='All genes')
    
    universal = comparison_df[comparison_df['target_type'] == 'universal']
    primary_spec = comparison_df[comparison_df['target_type'] == 'primary_specific']
    met_spec = comparison_df[comparison_df['target_type'] == 'metastatic_specific']
    
    if len(universal) > 0:
        ax.scatter(universal['exclusivity_mean_primary'],
                  universal['exclusivity_mean_metastatic'],
                  s=100, c='green', edgecolors='black', linewidth=1,
                  label='Universal', zorder=5)
        
        for _, row in universal.head(5).iterrows():
            ax.annotate(row['gene'],
                       (row['exclusivity_mean_primary'], row['exclusivity_mean_metastatic']),
                       fontsize=8, ha='left', va='bottom')
    
    if len(primary_spec) > 0:
        ax.scatter(primary_spec['exclusivity_mean_primary'],
                  primary_spec['exclusivity_mean_metastatic'],
                  s=80, c='blue', edgecolors='black', linewidth=1,
                  label='Primary-specific', zorder=4, alpha=0.7)
    
    if len(met_spec) > 0:
        ax.scatter(met_spec['exclusivity_mean_primary'],
                  met_spec['exclusivity_mean_metastatic'],
                  s=80, c='red', edgecolors='black', linewidth=1,
                  label='Metastatic-specific', zorder=4, alpha=0.7)
    
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


def create_subtype_comparison_heatmap(data, output_dir):
    """Figure 5: Subtype comparison heatmap"""
    print("\nCreating Figure 5: Subtype comparison...")
    
    subtype_comp = data['subtype_comparison']
    if subtype_comp is None:
        print("  No subtype comparison data - skipping")
        return
    
    available_subtypes = data['available_subtypes']
    if not available_subtypes:
        return
    
    # Get top 30 genes by mean exclusivity
    top_genes = subtype_comp.nlargest(30, 'mean_exclusivity')
    
    # Create matrix with only available subtypes
    matrix_data = []
    gene_names = []
    
    for _, row in top_genes.iterrows():
        gene_names.append(row['gene'])
        gene_values = [row.get(subtype, 0) for subtype in available_subtypes]
        matrix_data.append(gene_values)
    
    matrix = pd.DataFrame(matrix_data, columns=available_subtypes, index=gene_names)
    
    fig, ax = plt.subplots(figsize=(8, 12))
    
    sns.heatmap(matrix, cmap='RdYlBu_r', center=1.0, vmin=0, vmax=3,
               cbar_kws={'label': 'Stromal Exclusivity'},
               ax=ax, linewidths=0.5, linecolor='white',
               annot=False)
    
    ax.set_xlabel('Cancer Subtype', fontweight='bold')
    ax.set_ylabel('Gene', fontweight='bold')
    ax.set_title('Stromal Target Specificity Across Cancer Subtypes', 
                fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig5_subtype_comparison_heatmap.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig5_subtype_comparison_heatmap.pdf")


def create_subtype_barplot_comparison(data, output_dir):
    """Figure 6: Top targets per subtype (side-by-side bars)"""
    print("\nCreating Figure 6: Top targets per subtype...")
    
    available_subtypes = data['available_subtypes']
    if not available_subtypes:
        print("  No subtype data - skipping")
        return
    
    # Get top 10 genes per subtype
    top_genes_per_subtype = {}
    for subtype in available_subtypes:
        if data[subtype] is not None:
            top_10 = data[subtype].head(10)
            top_genes_per_subtype[subtype] = top_10[['gene', 'stromal_exclusivity']]
    
    if not top_genes_per_subtype:
        return
    
    # Create subplots
    n_subtypes = len(top_genes_per_subtype)
    fig, axes = plt.subplots(1, n_subtypes, figsize=(6*n_subtypes, 8), sharey=True)
    
    if n_subtypes == 1:
        axes = [axes]
    
    colors = {'ER+': 'pink', 'TNBC': 'purple', 'HER2+': 'orange', 'Normal': 'lightgreen'}
    
    for ax, (subtype, df) in zip(axes, top_genes_per_subtype.items()):
        y_pos = np.arange(len(df))
        ax.barh(y_pos, df['stromal_exclusivity'], color=colors.get(subtype, 'steelblue'))
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df['gene'])
        ax.invert_yaxis()
        ax.set_xlabel('Exclusivity', fontweight='bold')
        ax.set_title(f'{subtype} Top 10', fontweight='bold')
        ax.grid(alpha=0.3, axis='x')
    
    fig.suptitle('Top Stromal Targets by Cancer Subtype', 
                fontweight='bold', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig6_subtype_top_targets.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig6_subtype_top_targets.pdf")


def create_venn_subtype_overlap(data, output_dir):
    """Figure 7: Venn diagram showing subtype overlap"""
    print("\nCreating Figure 7: Subtype overlap...")
    
    subtype_comp = data['subtype_comparison']
    if subtype_comp is None:
        print("  No subtype comparison data - skipping")
        return
    
    available_subtypes = data['available_subtypes']
    if len(available_subtypes) < 2:
        print("  Need at least 2 subtypes for overlap analysis")
        return
    
    # Define high exclusivity threshold
    threshold = 0.5
    
    # Get sets of genes high in each subtype
    subtype_sets = {}
    for subtype in available_subtypes:
        if subtype in subtype_comp.columns:
            high_genes = subtype_comp[subtype_comp[subtype] > threshold]['gene'].tolist()
            subtype_sets[subtype] = set(high_genes)
    
    # Calculate overlaps
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if len(subtype_sets) == 2:
        # 2-way Venn
        try:
            from matplotlib_venn import venn2
            sets = list(subtype_sets.values())
            labels = list(subtype_sets.keys())
            venn2(sets, set_labels=labels, ax=ax)
            ax.set_title('Stromal Target Overlap Between Subtypes', 
                        fontweight='bold', fontsize=12)
        except ImportError:
            print("  matplotlib-venn not installed, creating bar plot instead")
            create_overlap_barplot(subtype_sets, ax)
    
    elif len(subtype_sets) == 3:
        # 3-way Venn
        try:
            from matplotlib_venn import venn3
            sets = list(subtype_sets.values())
            labels = list(subtype_sets.keys())
            venn3(sets, set_labels=labels, ax=ax)
            ax.set_title('Stromal Target Overlap Across Subtypes', 
                        fontweight='bold', fontsize=12)
        except ImportError:
            print("  matplotlib-venn not installed, creating bar plot instead")
            create_overlap_barplot(subtype_sets, ax)
    else:
        # More than 3 - use bar plot
        create_overlap_barplot(subtype_sets, ax)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig7_subtype_overlap.pdf', bbox_inches='tight')
    plt.close()
    
    print("  Saved: fig7_subtype_overlap.pdf")


def create_overlap_barplot(subtype_sets, ax):
    """Create bar plot showing unique and shared genes"""
    # Calculate unique and shared
    all_genes = set.union(*subtype_sets.values())
    
    overlaps = {}
    for gene in all_genes:
        present_in = [s for s, genes in subtype_sets.items() if gene in genes]
        if len(present_in) == 1:
            key = f"{present_in[0]} only"
        else:
            key = "Shared: " + ", ".join(sorted(present_in))
        overlaps[key] = overlaps.get(key, 0) + 1
    
    # Plot
    categories = list(overlaps.keys())
    counts = list(overlaps.values())
    
    y_pos = np.arange(len(categories))
    ax.barh(y_pos, counts, color='steelblue')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories)
    ax.set_xlabel('Number of Genes', fontweight='bold')
    ax.set_title('Stromal Target Overlap Across Subtypes', fontweight='bold')
    ax.grid(alpha=0.3, axis='x')


def create_summary_table(data, output_dir):
    """Create summary tables"""
    print("\nCreating summary tables...")
    
    # Global targets table
    if data['global'] is not None:
        top_targets = data['global'].head(30)[['gene', 'stromal_exclusivity', 'log2fc_vs_all',
                                               'stromal_detection', 'tumor_penalty']].copy()
        
        top_targets['stromal_exclusivity'] = top_targets['stromal_exclusivity'].round(3)
        top_targets['log2fc_vs_all'] = top_targets['log2fc_vs_all'].round(2)
        top_targets['stromal_detection'] = (top_targets['stromal_detection'] * 100).round(1)
        top_targets['tumor_penalty'] = top_targets['tumor_penalty'].round(2)
        
        top_targets.columns = ['Gene', 'Exclusivity', 'Log2FC', 'Detection (%)', 'Tumor Penalty']
        
        top_targets.to_csv(output_dir / 'table1_top_global_targets.csv', index=False)
        print("  Saved: table1_top_global_targets.csv")
    
    # Subtype comparison table
    if data['subtype_comparison'] is not None:
        subtype_table = data['subtype_comparison'].copy()
        
        # Round numeric columns
        numeric_cols = subtype_table.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            subtype_table[col] = subtype_table[col].round(3)
        
        subtype_table.head(50).to_csv(output_dir / 'table2_subtype_comparison.csv', index=False)
        print("  Saved: table2_subtype_comparison.csv")


def main():
    print("\n" + "="*80)
    print("STEP 5: CREATE VISUALIZATIONS (Subtype-Aware)")
    print("="*80)
    
    data = load_all_results()
    
    output_dir = Path('./results/figures')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Core figures
    create_top_targets_barplot(data, output_dir)
    create_volcano_plot(data, output_dir)
    # create_patient_heatmap(data, output_dir)
    create_primary_vs_metastatic_plot(data, output_dir)
    
    # Subtype-specific figures
    if data['available_subtypes']:
        create_subtype_comparison_heatmap(data, output_dir)
        create_subtype_barplot_comparison(data, output_dir)
        create_venn_subtype_overlap(data, output_dir)
    
    # Summary tables
    create_summary_table(data, output_dir)
    
    print("\n" + "="*80)
    print("STEP 5 COMPLETE")
    print("="*80)
    print(f"\nAll figures saved in: {output_dir.absolute()}")
    
    print("\nFigures created:")
    print("  1. fig1_top_global_targets.pdf       - Top 20 global targets")
    print("  2. fig2_volcano.pdf                   - Specificity landscape")
    print("  3. fig3_patient_heatmap.pdf           - Patient heterogeneity")
    print("  4. fig4_primary_vs_metastatic.pdf     - Stage comparison")
    
    if data['available_subtypes']:
        print("  5. fig5_subtype_comparison_heatmap.pdf - Subtype comparison")
        print("  6. fig6_subtype_top_targets.pdf        - Top targets per subtype")
        print("  7. fig7_subtype_overlap.pdf            - Subtype overlap analysis")
    
    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
