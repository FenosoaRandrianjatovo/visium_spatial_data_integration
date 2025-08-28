#!/usr/bin/env python3
"""
Spatial Data Integration Pipeline using SpatialGlue
=================================================

This script performs spatial data integration using SpatialGlue on the 
preprocessed Visium data.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Import SpatialGlue
SPATIALGLUE_AVAILABLE = False
try:
    import spatialglue as sg
    SPATIALGLUE_AVAILABLE = True
    print("✅ SpatialGlue imported successfully")
except ImportError as e:
    print(f"⚠️ SpatialGlue not available: {e}")
    print("Will use alternative integration methods")

# Set settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_processed_data():
    """Load the preprocessed data files."""
    print("🔄 Loading processed data...")
    
    data_files = {
        'spatial': 'processed_data/visium_spatial_processed.h5ad',
        'spatial_abs': 'processed_data/visium_with_abs_processed.h5ad'
    }
    
    loaded_data = {}
    
    for name, filepath in data_files.items():
        if os.path.exists(filepath):
            adata = sc.read_h5ad(filepath)
            loaded_data[name] = adata
            print(f"✅ Loaded {name}: {adata.shape}")
        else:
            print(f"⚠️ File not found: {filepath}")
    
    return loaded_data

def prepare_for_integration(adata_dict):
    """
    Prepare data for SpatialGlue integration.
    
    Args:
        adata_dict: Dictionary of AnnData objects
        
    Returns:
        List of prepared AnnData objects
    """
    print("🔄 Preparing data for integration...")
    
    prepared_data = []
    
    for name, adata in adata_dict.items():
        print(f"📊 Preparing {name}...")
        
        # Copy to avoid modifying original
        adata_prep = adata.copy()
        
        # Ensure we have highly variable genes
        if 'highly_variable' not in adata_prep.var.columns:
            sc.pp.highly_variable_genes(adata_prep, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        # Filter to highly variable genes for integration
        adata_prep = adata_prep[:, adata_prep.var.highly_variable].copy()
        
        # Scale data
        sc.pp.scale(adata_prep, max_value=10)
        
        # Add slice name
        adata_prep.obs['slice_name'] = name
        
        prepared_data.append(adata_prep)
        print(f"✅ Prepared {name}: {adata_prep.shape}")
    
    return prepared_data

def run_spatialglue_integration(adata_list):
    """
    Run SpatialGlue integration on the prepared data.
    
    Args:
        adata_list: List of AnnData objects
        
    Returns:
        Integrated AnnData object or None if SpatialGlue is not available
    """
    if not SPATIALGLUE_AVAILABLE:
        print("⚠️ SpatialGlue not available, skipping...")
        return None
        
    print("🔄 Running SpatialGlue integration...")
    
    try:
        # Create slice names
        slice_names = [adata.obs['slice_name'].iloc[0] for adata in adata_list]
        print(f"📊 Integrating slices: {slice_names}")
        
        # Run SpatialGlue integration
        # Note: SpatialGlue parameters may need adjustment based on your specific data
        integrated_adata = sg.spatial_integration(
            adata_list,
            slice_names=slice_names,
            # Common parameters - adjust as needed
            n_neighbors=15,
            n_pcs=50,
            resolution=0.5,
            # You may need to adjust these based on SpatialGlue documentation
        )
        
        print(f"✅ Integration completed: {integrated_adata.shape}")
        return integrated_adata
        
    except Exception as e:
        print(f"❌ Error during SpatialGlue integration: {e}")
        print("Falling back to alternative method...")
        return None

def alternative_integration_approach(adata_list):
    """
    Alternative integration approach using scanpy/scvi-tools methods.
    This serves as a backup if SpatialGlue has issues.
    """
    print("🔄 Running alternative integration with scanpy...")
    
    # Concatenate the datasets
    adata_concat = sc.concat(adata_list, join='outer', index_unique='_')
    
    # Add batch information
    adata_concat.obs['batch'] = adata_concat.obs['slice_name']
    
    # Handle NaN values
    print("🔄 Handling missing values...")
    import numpy as np
    
    # Replace NaN in X matrix with 0
    if np.isnan(adata_concat.X).any():
        print("⚠️ Found NaN values in expression matrix, replacing with 0")
        adata_concat.X = np.nan_to_num(adata_concat.X, nan=0.0)
    
    # Handle NaN in obs columns
    for col in adata_concat.obs.columns:
        if adata_concat.obs[col].isna().any():
            if adata_concat.obs[col].dtype.name == 'category':
                # Handle categorical columns
                adata_concat.obs[col] = adata_concat.obs[col].cat.add_categories(['Unknown'])
                adata_concat.obs[col] = adata_concat.obs[col].fillna('Unknown')
            elif adata_concat.obs[col].dtype == 'object':
                adata_concat.obs[col] = adata_concat.obs[col].fillna('Unknown')
            else:
                adata_concat.obs[col] = adata_concat.obs[col].fillna(0)
    
    # Find highly variable genes across batches
    sc.pp.highly_variable_genes(adata_concat, batch_key='batch', n_top_genes=2000)
    
    # Keep only highly variable genes
    adata_concat = adata_concat[:, adata_concat.var.highly_variable]
    
    # Scale the data
    sc.pp.scale(adata_concat, max_value=10)
    
    # PCA
    sc.tl.pca(adata_concat, svd_solver='arpack', n_comps=50)
    
    # Batch correction using Harmony (if available) or basic integration
    try:
        import harmonypy as hm
        harmony_out = hm.run_harmony(
            adata_concat.obsm['X_pca'], 
            adata_concat.obs, 
            vars_use=['batch']
        )
        adata_concat.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
        
        # Compute neighborhood graph
        sc.pp.neighbors(adata_concat, use_rep='X_pca_harmony', n_neighbors=15)
        print("✅ Used Harmony for batch correction")
        
    except ImportError:
        # Fallback to basic integration
        sc.pp.neighbors(adata_concat, n_neighbors=15)
        print("⚠️ Harmony not available, using basic integration")
    
    # UMAP
    sc.tl.umap(adata_concat)
    
    # Leiden clustering
    sc.tl.leiden(adata_concat, resolution=0.5)
    
    print(f"✅ Alternative integration completed: {adata_concat.shape}")
    return adata_concat

def create_integration_plots(adata, output_prefix):
    """Create plots to visualize the integration results."""
    print(f"📊 Creating integration plots...")
    
    # Create output directory for integration plots
    os.makedirs('plots/integration', exist_ok=True)
    
    # Plot 1: UMAP colored by batch/slice
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Spatial Data Integration Results - {output_prefix}', fontsize=16)
    
    # UMAP by batch
    if 'slice_name' in adata.obs.columns:
        sc.pl.umap(adata, color='slice_name', ax=axes[0,0], show=False, frameon=False)
        axes[0,0].set_title('Integration by Dataset')
    elif 'batch' in adata.obs.columns:
        sc.pl.umap(adata, color='batch', ax=axes[0,0], show=False, frameon=False)
        axes[0,0].set_title('Integration by Batch')
    
    # UMAP by zonation
    if 'zonationGroup' in adata.obs.columns:
        sc.pl.umap(adata, color='zonationGroup', ax=axes[0,1], show=False, frameon=False)
        axes[0,1].set_title('Zonation Groups')
    
    # UMAP by clusters
    if 'leiden' in adata.obs.columns:
        sc.pl.umap(adata, color='leiden', ax=axes[0,2], show=False, frameon=False)
        axes[0,2].set_title('Leiden Clusters')
    
    # Sample distribution
    if 'sample' in adata.obs.columns:
        adata.obs['sample'].value_counts().plot(kind='bar', ax=axes[1,0])
        axes[1,0].set_title('Spots per Sample')
        axes[1,0].tick_params(axis='x', rotation=45)
    elif 'orig.ident' in adata.obs.columns:
        adata.obs['orig.ident'].value_counts().plot(kind='bar', ax=axes[1,0])
        axes[1,0].set_title('Spots per Sample')
        axes[1,0].tick_params(axis='x', rotation=45)
    
    # Dataset distribution
    if 'slice_name' in adata.obs.columns:
        adata.obs['slice_name'].value_counts().plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Spots per Dataset')
        axes[1,1].tick_params(axis='x', rotation=45)
    elif 'batch' in adata.obs.columns:
        adata.obs['batch'].value_counts().plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Spots per Batch')
        axes[1,1].tick_params(axis='x', rotation=45)
    
    # Cluster distribution
    if 'leiden' in adata.obs.columns:
        adata.obs['leiden'].value_counts().plot(kind='bar', ax=axes[1,2])
        axes[1,2].set_title('Spots per Cluster')
        axes[1,2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f'plots/integration/{output_prefix}_integration_overview.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create separate detailed UMAP plots
    if 'leiden' in adata.obs.columns:
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='leiden', legend_loc='on data', 
                   legend_fontsize=8, show=False, frameon=False)
        plt.title('Integrated Clusters (Leiden)')
        plt.savefig(f'plots/integration/{output_prefix}_clusters_detailed.png', 
                    dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"✅ Integration plots saved")

def save_integration_results(adata, filename):
    """Save integration results."""
    os.makedirs('integrated_data', exist_ok=True)
    filepath = f'integrated_data/{filename}'
    
    # Remove raw data to save space
    adata_copy = adata.copy()
    adata_copy.raw = None
    
    adata_copy.write(filepath)
    print(f"✅ Saved integrated data: {filepath}")

def main():
    """Main integration pipeline."""
    print("🚀 Starting Spatial Data Integration Pipeline")
    print("=" * 50)
    
    # Load processed data
    data_dict = load_processed_data()
    
    if len(data_dict) < 2:
        print("❌ Need at least 2 datasets for integration")
        print("Please run data_processing.py first")
        return
    
    # Prepare data for integration
    prepared_data = prepare_for_integration(data_dict)
    
    # Try SpatialGlue integration first
    integrated_adata = run_spatialglue_integration(prepared_data)
    
    # If SpatialGlue fails, use alternative approach
    if integrated_adata is None:
        print("🔄 Falling back to alternative integration method...")
        integrated_adata = alternative_integration_approach(prepared_data)
    
    if integrated_adata is not None:
        # Create plots
        create_integration_plots(integrated_adata, 'spatialglue')
        
        # Save results
        save_integration_results(integrated_adata, 'integrated_visium_data.h5ad')
        
        print("\n🎉 Spatial integration completed!")
        print("📁 Check the following directories:")
        print("   - integrated_data/ : Integrated H5AD file")
        print("   - plots/integration/ : Integration visualization plots")
        
        # Print summary
        print(f"\n📊 Integration Summary:")
        print(f"   - Final shape: {integrated_adata.shape}")
        print(f"   - Datasets integrated: {integrated_adata.obs['slice_name'].nunique() if 'slice_name' in integrated_adata.obs else integrated_adata.obs['batch'].nunique()}")
        print(f"   - Total clusters: {integrated_adata.obs['leiden'].nunique() if 'leiden' in integrated_adata.obs else 'N/A'}")
        
    else:
        print("❌ Integration failed with both methods")

if __name__ == "__main__":
    main()
