#!/usr/bin/env python3
"""
Quick Integration Test
=====================
A simplified version to test integration without saving large files.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

def main():
    print("🔄 Quick integration test...")
    
    # Load processed data
    adata_spatial = sc.read_h5ad('processed_data/visium_spatial_processed.h5ad')
    adata_abs = sc.read_h5ad('processed_data/visium_with_abs_processed.h5ad')
    
    print(f"Spatial: {adata_spatial.shape}")
    print(f"With ABs: {adata_abs.shape}")
    
    # Quick prep
    adata_spatial.obs['slice_name'] = 'spatial'
    adata_abs.obs['slice_name'] = 'spatial_abs'
    
    # Filter to top variable genes only
    n_genes = min(500, adata_spatial.var.highly_variable.sum(), adata_abs.var.highly_variable.sum())
    
    adata_spatial_subset = adata_spatial[:, adata_spatial.var.highly_variable][:, :n_genes].copy()
    adata_abs_subset = adata_abs[:, adata_abs.var.highly_variable][:, :n_genes].copy()
    
    # Concatenate
    adata_concat = sc.concat([adata_spatial_subset, adata_abs_subset], join='outer', index_unique='_')
    adata_concat.obs['batch'] = adata_concat.obs['slice_name']
    
    # Handle missing values
    adata_concat.X = np.nan_to_num(adata_concat.X, nan=0.0)
    
    # Basic integration
    sc.pp.scale(adata_concat, max_value=10)
    sc.tl.pca(adata_concat, n_comps=30)
    sc.pp.neighbors(adata_concat, n_neighbors=15)
    sc.tl.umap(adata_concat)
    sc.tl.leiden(adata_concat, resolution=0.5)
    
    print(f"✅ Integration completed: {adata_concat.shape}")
    print(f"Clusters: {adata_concat.obs['leiden'].nunique()}")
    
    # Create a simple plot
    os.makedirs('plots/integration', exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot by batch
    for i, batch in enumerate(adata_concat.obs['batch'].unique()):
        mask = adata_concat.obs['batch'] == batch
        axes[0].scatter(adata_concat.obsm['X_umap'][mask, 0], 
                       adata_concat.obsm['X_umap'][mask, 1],
                       label=batch, alpha=0.7, s=1)
    axes[0].legend()
    axes[0].set_title('Integration by Dataset')
    
    # Plot by clusters
    clusters = adata_concat.obs['leiden'].astype('category').cat.codes
    scatter = axes[1].scatter(adata_concat.obsm['X_umap'][:, 0], 
                             adata_concat.obsm['X_umap'][:, 1],
                             c=clusters, cmap='tab20', alpha=0.7, s=1)
    axes[1].set_title('Leiden Clusters')
    
    plt.tight_layout()
    plt.savefig('plots/integration/quick_integration_test.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("✅ Integration test completed!")
    print("📊 Summary:")
    print(f"   - Total cells: {adata_concat.shape[0]}")
    print(f"   - Datasets: {adata_concat.obs['batch'].nunique()}")
    print(f"   - Clusters: {adata_concat.obs['leiden'].nunique()}")
    print(f"   - Plot saved: plots/integration/quick_integration_test.png")

if __name__ == "__main__":
    main()
