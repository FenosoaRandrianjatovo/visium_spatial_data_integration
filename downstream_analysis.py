#!/usr/bin/env python3
"""
Downstream Analysis Pipeline
===========================

This script performs downstream analysis on the integrated spatial data,
including differential expression, pathway analysis, and spatial patterns.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Set settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_integrated_data():
    """Load the integrated data."""
    print("🔄 Loading integrated data...")
    
    filepath = 'integrated_data/integrated_visium_data.h5ad'
    if os.path.exists(filepath):
        adata = sc.read_h5ad(filepath)
        print(f"✅ Loaded integrated data: {adata.shape}")
        return adata
    else:
        print(f"❌ Integrated data file not found: {filepath}")
        print("Please run spatial_integration.py first")
        return None

def differential_expression_analysis(adata):
    """Perform differential expression analysis."""
    print("🔄 Running differential expression analysis...")
    
    # Create output directory
    os.makedirs('results/differential_expression', exist_ok=True)
    
    # DE analysis by clusters
    if 'leiden' in adata.obs.columns:
        print("📊 DE analysis by clusters...")
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
        
        # Save results
        de_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        de_df.to_csv('results/differential_expression/cluster_markers.csv', index=False)
        
        # Plot top markers
        sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, show=False)
        plt.savefig('results/differential_expression/cluster_markers_plot.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Cluster DE analysis completed")
    
    # DE analysis by zonation groups
    if 'zonationGroup' in adata.obs.columns:
        print("📊 DE analysis by zonation...")
        sc.tl.rank_genes_groups(adata, 'zonationGroup', method='wilcoxon')
        
        # Save results
        de_zonation_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        de_zonation_df.to_csv('results/differential_expression/zonation_markers.csv', index=False)
        
        # Plot top markers
        sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, show=False)
        plt.savefig('results/differential_expression/zonation_markers_plot.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Zonation DE analysis completed")
    
    # DE analysis by dataset
    if 'slice_name' in adata.obs.columns:
        print("📊 DE analysis by dataset...")
        sc.tl.rank_genes_groups(adata, 'slice_name', method='wilcoxon')
        
        # Save results
        de_dataset_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        de_dataset_df.to_csv('results/differential_expression/dataset_markers.csv', index=False)
        
        print("✅ Dataset DE analysis completed")

def analyze_spatial_patterns(adata):
    """Analyze spatial patterns in the data."""
    print("🔄 Analyzing spatial patterns...")
    
    # Create output directory
    os.makedirs('results/spatial_analysis', exist_ok=True)
    
    # If we have original spatial coordinates
    if 'UMAP_1' in adata.obs.columns and 'UMAP_2' in adata.obs.columns:
        print("📊 Creating spatial pattern plots...")
        
        # Plot spatial distribution of clusters
        if 'leiden' in adata.obs.columns:
            plt.figure(figsize=(12, 8))
            scatter = plt.scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'], 
                                c=adata.obs['leiden'].astype('category').cat.codes,
                                cmap='tab20', alpha=0.7, s=2)
            plt.colorbar(scatter)
            plt.title('Spatial Distribution of Clusters')
            plt.xlabel('UMAP_1')
            plt.ylabel('UMAP_2')
            plt.savefig('results/spatial_analysis/spatial_clusters.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Plot spatial distribution of zonation
        if 'zonationGroup' in adata.obs.columns:
            plt.figure(figsize=(12, 8))
            zones = adata.obs['zonationGroup'].unique()
            colors = plt.cm.viridis(np.linspace(0, 1, len(zones)))
            
            for i, zone in enumerate(zones):
                mask = adata.obs['zonationGroup'] == zone
                plt.scatter(adata.obs.loc[mask, 'UMAP_1'], 
                           adata.obs.loc[mask, 'UMAP_2'],
                           c=[colors[i]], label=zone, alpha=0.7, s=2)
            
            plt.legend()
            plt.title('Spatial Distribution of Zonation Groups')
            plt.xlabel('UMAP_1')
            plt.ylabel('UMAP_2')
            plt.savefig('results/spatial_analysis/spatial_zonation.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    # Gene expression spatial patterns
    print("📊 Analyzing gene expression patterns...")
    
    # Select top variable genes for spatial analysis
    if hasattr(adata, 'raw') and adata.raw is not None:
        top_genes = adata.var.nlargest(20, 'dispersions_norm').index
    else:
        # Fallback: use genes with highest mean expression
        gene_means = np.array(adata.X.mean(axis=0)).flatten()
        top_gene_indices = np.argsort(gene_means)[-20:]
        top_genes = adata.var.index[top_gene_indices]
    
    # Plot spatial expression of top genes
    if 'UMAP_1' in adata.obs.columns and 'UMAP_2' in adata.obs.columns:
        n_genes_plot = min(9, len(top_genes))
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        axes = axes.flatten()
        
        for i, gene in enumerate(top_genes[:n_genes_plot]):
            if gene in adata.var.index:
                gene_expr = adata[:, gene].X
                if hasattr(gene_expr, 'toarray'):
                    gene_expr = gene_expr.toarray().flatten()
                else:
                    gene_expr = gene_expr.flatten()
                
                scatter = axes[i].scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'],
                                        c=gene_expr, cmap='viridis', alpha=0.7, s=1)
                axes[i].set_title(f'{gene}')
                axes[i].set_xlabel('UMAP_1')
                axes[i].set_ylabel('UMAP_2')
                plt.colorbar(scatter, ax=axes[i])
        
        plt.tight_layout()
        plt.savefig('results/spatial_analysis/top_genes_spatial.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    print("✅ Spatial pattern analysis completed")

def create_summary_report(adata):
    """Create a comprehensive summary report."""
    print("📊 Creating summary report...")
    
    # Create output directory
    os.makedirs('results/summary', exist_ok=True)
    
    # Collect summary statistics
    summary_stats = {
        'Total spots': adata.shape[0],
        'Total genes': adata.shape[1],
        'Datasets': adata.obs['slice_name'].nunique() if 'slice_name' in adata.obs else adata.obs['batch'].nunique() if 'batch' in adata.obs else 'N/A',
        'Clusters': adata.obs['leiden'].nunique() if 'leiden' in adata.obs else 'N/A',
        'Zonation groups': adata.obs['zonationGroup'].nunique() if 'zonationGroup' in adata.obs else 'N/A',
        'Samples': adata.obs['sample'].nunique() if 'sample' in adata.obs else adata.obs['orig.ident'].nunique() if 'orig.ident' in adata.obs else 'N/A'
    }
    
    # Create summary plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Integrated Spatial Data Summary', fontsize=16)
    
    # Plot 1: Data overview
    categories = list(summary_stats.keys())
    values = [v if isinstance(v, (int, float)) else 0 for v in summary_stats.values()]
    axes[0,0].bar(categories[:4], values[:4])
    axes[0,0].set_title('Data Overview')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Cluster distribution
    if 'leiden' in adata.obs.columns:
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        axes[0,1].bar(cluster_counts.index, cluster_counts.values)
        axes[0,1].set_title('Cluster Distribution')
        axes[0,1].set_xlabel('Cluster')
        axes[0,1].set_ylabel('Number of spots')
    
    # Plot 3: Dataset distribution
    if 'slice_name' in adata.obs.columns:
        dataset_counts = adata.obs['slice_name'].value_counts()
        axes[0,2].pie(dataset_counts.values, labels=dataset_counts.index, autopct='%1.1f%%')
        axes[0,2].set_title('Dataset Distribution')
    elif 'batch' in adata.obs.columns:
        batch_counts = adata.obs['batch'].value_counts()
        axes[0,2].pie(batch_counts.values, labels=batch_counts.index, autopct='%1.1f%%')
        axes[0,2].set_title('Batch Distribution')
    
    # Plot 4: Zonation distribution
    if 'zonationGroup' in adata.obs.columns:
        zone_counts = adata.obs['zonationGroup'].value_counts()
        axes[1,0].bar(zone_counts.index, zone_counts.values)
        axes[1,0].set_title('Zonation Distribution')
        axes[1,0].tick_params(axis='x', rotation=45)
    
    # Plot 5: Quality metrics
    if 'total_counts' in adata.obs.columns:
        axes[1,1].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
        axes[1,1].set_title('Total Counts per Spot')
        axes[1,1].set_xlabel('Total counts')
        axes[1,1].set_ylabel('Number of spots')
    
    # Plot 6: Gene counts
    if 'n_genes_by_counts' in adata.obs.columns:
        axes[1,2].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7)
        axes[1,2].set_title('Genes per Spot')
        axes[1,2].set_xlabel('Number of genes')
        axes[1,2].set_ylabel('Number of spots')
    
    plt.tight_layout()
    plt.savefig('results/summary/integration_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save summary statistics
    summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Value'])
    summary_df.to_csv('results/summary/summary_statistics.csv')
    
    # Create detailed summary text
    with open('results/summary/analysis_summary.txt', 'w') as f:
        f.write("SPATIAL DATA INTEGRATION ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("DATASET OVERVIEW:\n")
        f.write("-" * 20 + "\n")
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\n")
        
        f.write("\nFILES GENERATED:\n")
        f.write("-" * 20 + "\n")
        f.write("1. processed_data/ - Preprocessed individual datasets\n")
        f.write("2. integrated_data/ - Integrated dataset\n")
        f.write("3. plots/ - Quality control and integration plots\n")
        f.write("4. results/differential_expression/ - DE analysis results\n")
        f.write("5. results/spatial_analysis/ - Spatial pattern analysis\n")
        f.write("6. results/summary/ - Summary statistics and plots\n")
        
        if 'leiden' in adata.obs.columns:
            f.write(f"\nCLUSTER INFORMATION:\n")
            f.write("-" * 20 + "\n")
            cluster_counts = adata.obs['leiden'].value_counts().sort_index()
            for cluster, count in cluster_counts.items():
                f.write(f"Cluster {cluster}: {count} spots\n")
        
        if 'zonationGroup' in adata.obs.columns:
            f.write(f"\nZONATION INFORMATION:\n")
            f.write("-" * 20 + "\n")
            zone_counts = adata.obs['zonationGroup'].value_counts()
            for zone, count in zone_counts.items():
                f.write(f"{zone}: {count} spots\n")
    
    print("✅ Summary report created")

def main():
    """Main downstream analysis pipeline."""
    print("🚀 Starting Downstream Analysis Pipeline")
    print("=" * 50)
    
    # Load integrated data
    adata = load_integrated_data()
    if adata is None:
        return
    
    # Create results directory
    os.makedirs('results', exist_ok=True)
    
    # Run differential expression analysis
    differential_expression_analysis(adata)
    
    # Analyze spatial patterns
    analyze_spatial_patterns(adata)
    
    # Create summary report
    create_summary_report(adata)
    
    print("\n🎉 Downstream analysis completed!")
    print("📁 Check the 'results/' directory for all analysis outputs")

if __name__ == "__main__":
    main()
