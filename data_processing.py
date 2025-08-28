#!/usr/bin/env python3
"""
Visium Spatial Data Processing Pipeline
=====================================

This script processes Visium spatial transcriptomics data and prepares it for 
spatial integration using SpatialGlue.

Requirements:
- scanpy
- pandas 
- anndata
- spatialglue
- matplotlib
- seaborn
"""

import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from pathlib import Path

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def setup_directories():
    """Create output directories for processed data and plots."""
    directories = ['processed_data', 'plots', 'logs']
    for dir_name in directories:
        os.makedirs(dir_name, exist_ok=True)
    print("✅ Output directories created")

def load_visium_spatial_data():
    """
    Load and process the standard Visium spatial data.
    
    Returns:
        AnnData: Processed spatial transcriptomics data
    """
    print("\n🔄 Loading Visium Spatial Data...")
    
    # Load the 10x formatted data
    try:
        # Read the files manually since the feature file has non-standard format
        import scipy.io
        import gzip
        
        # Read the matrix manually
        matrix = scipy.io.mmread('Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz')
        
        # Read features manually
        with gzip.open('Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/features.tsv.gz', 'rt') as f:
            features = [line.strip() for line in f]
        
        # Read barcodes manually
        with gzip.open('Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/barcodes.tsv.gz', 'rt') as f:
            barcodes = [line.strip() for line in f]
        
        # Create AnnData object
        adata = ad.AnnData(
            X=matrix.T.tocsr(),  # Transpose and convert to CSR format (cells x genes)
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=features)
        )
        
        adata.var_names_make_unique()
        print(f"✅ Loaded count matrix: {adata.shape}")
        
    except Exception as e:
        print(f"❌ Error loading count matrix: {e}")
        return None
    
    # Load annotations
    try:
        annotations = pd.read_csv('Data/Visium_Spatial/annot_mouseStStVisium.csv')
        print(f"✅ Loaded annotations: {annotations.shape}")
        
        # Set barcode as index for merging
        annotations = annotations.set_index('spot')
        
        # Add annotations to adata
        # Filter adata to only include spots in annotations
        common_spots = adata.obs.index.intersection(annotations.index)
        adata = adata[common_spots].copy()
        
        # Add annotation data
        for col in annotations.columns:
            adata.obs[col] = annotations.loc[adata.obs.index, col]
            
        print(f"✅ Merged annotations. Final shape: {adata.shape}")
        
    except Exception as e:
        print(f"❌ Error loading annotations: {e}")
        return None
    
    # Add metadata
    adata.obs['dataset'] = 'visium_spatial'
    adata.obs['technology'] = 'visium'
    
    return adata

def load_visium_with_abs_data():
    """
    Load and process the Visium spatial data with antibodies.
    
    Returns:
        AnnData: Processed spatial transcriptomics + protein data
    """
    print("\n🔄 Loading Visium + Antibodies Data...")
    
    try:
        # Load the h5 file using scanpy's 10x h5 reader
        adata = sc.read_10x_h5('Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_mouseStStVisiumWithABs.h5')
        adata.var_names_make_unique()
        # Note: Don't transpose for 10x h5 format - cells are already rows
        print(f"✅ Loaded h5 data: {adata.shape}")
        
        # Load annotations
        annotations = pd.read_csv('Data/Visium_spatial_with_ABs/annot_mouseStStWithABsVisium.csv')
        annotations = annotations.set_index('spot')
        print(f"✅ Loaded annotations: {annotations.shape}")
        
        # Filter and merge
        common_spots = adata.obs.index.intersection(annotations.index)
        print(f"🔍 Common spots found: {len(common_spots)}")
        
        if len(common_spots) == 0:
            print("⚠️ No common spots found. Checking barcode format...")
            # Try to match without the sample suffix
            adata_barcodes_clean = [bc.split('-')[0] for bc in adata.obs.index]
            annot_barcodes_clean = [bc.split('-')[0] for bc in annotations.index]
            
            # Find intersection with clean barcodes
            clean_common = set(adata_barcodes_clean).intersection(set(annot_barcodes_clean))
            print(f"🔍 Common spots with clean barcodes: {len(clean_common)}")
            
            if len(clean_common) > 0:
                # Create mapping
                adata_mapping = {bc.split('-')[0]: bc for bc in adata.obs.index}
                annot_mapping = {bc.split('-')[0]: bc for bc in annotations.index}
                
                # Filter to common clean barcodes
                common_clean_list = list(clean_common)
                adata_subset = [adata_mapping[bc] for bc in common_clean_list]
                annot_subset = [annot_mapping[bc] for bc in common_clean_list]
                
                # Subset the data
                adata = adata[adata_subset].copy()
                annotations_subset = annotations.loc[annot_subset].copy()
                annotations_subset.index = adata.obs.index  # Match indices
                
                # Add annotation data
                for col in annotations_subset.columns:
                    adata.obs[col] = annotations_subset[col]
                    
                print(f"✅ Merged with clean barcodes. Final shape: {adata.shape}")
            else:
                print("❌ No matching barcodes found even with cleaned format")
                return None
        else:
            adata = adata[common_spots].copy()
            
            # Add annotation data
            for col in annotations.columns:
                adata.obs[col] = annotations.loc[adata.obs.index, col]
                
            print(f"✅ Merged annotations. Final shape: {adata.shape}")
        
        # Load ADT data separately
        try:
            adt_data = pd.read_csv('Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_ADT_mouseStStVisiumWithABs.csv', 
                                 index_col=0)
            print(f"✅ Loaded ADT data: {adt_data.shape}")
            
            # Add ADT data as a separate layer or obsm
            common_adt_spots = adata.obs.index.intersection(adt_data.index)
            if len(common_adt_spots) > 0:
                adata.obsm['protein'] = adt_data.loc[common_adt_spots]
                print(f"✅ Added protein data for {len(common_adt_spots)} spots")
                
        except Exception as e:
            print(f"⚠️ Could not load ADT data: {e}")
        
    except Exception as e:
        print(f"❌ Error loading Visium+ABs data: {e}")
        return None
    
    # Add metadata
    adata.obs['dataset'] = 'visium_with_abs'
    adata.obs['technology'] = 'visium_protein'
    
    return adata

def basic_preprocessing(adata, min_genes=200, min_cells=3, target_sum=1e4):
    """
    Perform basic preprocessing on the data.
    
    Args:
        adata: AnnData object
        min_genes: Minimum number of genes per cell
        min_cells: Minimum number of cells per gene
        target_sum: Target sum for normalization
        
    Returns:
        AnnData: Preprocessed data
    """
    if adata.shape[0] == 0:
        print("❌ Empty dataset - cannot preprocess")
        return adata
        
    dataset_name = adata.obs['dataset'].iloc[0] if 'dataset' in adata.obs.columns and len(adata.obs) > 0 else 'unknown'
    print(f"\n🔄 Preprocessing {dataset_name}...")
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    print(f"📊 Before filtering: {adata.shape}")
    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"📊 After filtering: {adata.shape}")
    
    # Save raw data
    adata.raw = adata
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    print(f"✅ Found {adata.var.highly_variable.sum()} highly variable genes")
    
    return adata

def create_qc_plots(adata, output_prefix):
    """Create quality control plots."""
    print(f"📊 Creating QC plots for {output_prefix}...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Quality Control Metrics - {output_prefix}', fontsize=16)
    
    # Plot 1: Number of genes per spot
    sns.histplot(adata.obs['n_genes_by_counts'], bins=50, ax=axes[0,0])
    axes[0,0].set_title('Genes per spot')
    axes[0,0].axvline(x=200, color='red', linestyle='--', alpha=0.7)
    
    # Plot 2: Total counts per spot
    sns.histplot(adata.obs['total_counts'], bins=50, ax=axes[0,1])
    axes[0,1].set_title('Total counts per spot')
    
    # Plot 3: Mitochondrial gene percentage
    if 'pct_counts_mt' in adata.obs.columns:
        sns.histplot(adata.obs['pct_counts_mt'], bins=50, ax=axes[0,2])
        axes[0,2].set_title('Mitochondrial gene %')
        axes[0,2].axvline(x=20, color='red', linestyle='--', alpha=0.7)
    
    # Plot 4: Sample distribution
    if 'sample' in adata.obs.columns:
        adata.obs['sample'].value_counts().plot(kind='bar', ax=axes[1,0])
        axes[1,0].set_title('Spots per sample')
        axes[1,0].tick_params(axis='x', rotation=45)
    elif 'orig.ident' in adata.obs.columns:
        adata.obs['orig.ident'].value_counts().plot(kind='bar', ax=axes[1,0])
        axes[1,0].set_title('Spots per sample')
        axes[1,0].tick_params(axis='x', rotation=45)
    
    # Plot 5: Zonation distribution
    if 'zonationGroup' in adata.obs.columns:
        adata.obs['zonationGroup'].value_counts().plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Spots per zone')
        axes[1,1].tick_params(axis='x', rotation=45)
    
    # Plot 6: Cluster distribution
    if 'cluster' in adata.obs.columns:
        adata.obs['cluster'].value_counts().plot(kind='bar', ax=axes[1,2])
        axes[1,2].set_title('Spots per cluster')
        axes[1,2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f'plots/{output_prefix}_qc_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create UMAP plot if coordinates exist
    if 'UMAP_1' in adata.obs.columns and 'UMAP_2' in adata.obs.columns:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        # UMAP colored by zonation
        if 'zonationGroup' in adata.obs.columns:
            scatter = axes[0].scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'], 
                                    c=adata.obs['zonationGroup'].astype('category').cat.codes, 
                                    cmap='viridis', alpha=0.7, s=1)
            axes[0].set_title('Zonation Groups')
            axes[0].set_xlabel('UMAP_1')
            axes[0].set_ylabel('UMAP_2')
            
        # UMAP colored by sample
        if 'sample' in adata.obs.columns:
            scatter = axes[1].scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'], 
                                    c=adata.obs['sample'].astype('category').cat.codes, 
                                    cmap='tab10', alpha=0.7, s=1)
            axes[1].set_title('Samples')
            axes[1].set_xlabel('UMAP_1')
            axes[1].set_ylabel('UMAP_2')
        elif 'orig.ident' in adata.obs.columns:
            scatter = axes[1].scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'], 
                                    c=adata.obs['orig.ident'].astype('category').cat.codes, 
                                    cmap='tab10', alpha=0.7, s=1)
            axes[1].set_title('Samples')
            axes[1].set_xlabel('UMAP_1')
            axes[1].set_ylabel('UMAP_2')
            
        # UMAP colored by cluster
        if 'cluster' in adata.obs.columns:
            scatter = axes[2].scatter(adata.obs['UMAP_1'], adata.obs['UMAP_2'], 
                                    c=adata.obs['cluster'].astype('category').cat.codes, 
                                    cmap='tab20', alpha=0.7, s=1)
            axes[2].set_title('Clusters')
            axes[2].set_xlabel('UMAP_1')
            axes[2].set_ylabel('UMAP_2')
            
        plt.tight_layout()
        plt.savefig(f'plots/{output_prefix}_umap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"✅ QC plots saved for {output_prefix}")

def save_processed_data(adata, filename):
    """Save processed data as h5ad file."""
    filepath = f'processed_data/{filename}'
    adata.write(filepath)
    print(f"✅ Saved processed data: {filepath}")

def main():
    """Main processing pipeline."""
    print("🚀 Starting Visium Data Processing Pipeline")
    print("=" * 50)
    
    # Setup
    setup_directories()
    
    # Process Visium Spatial data
    adata_spatial = load_visium_spatial_data()
    if adata_spatial is not None:
        adata_spatial = basic_preprocessing(adata_spatial)
        create_qc_plots(adata_spatial, 'visium_spatial')
        save_processed_data(adata_spatial, 'visium_spatial_processed.h5ad')
    
    # Process Visium + Antibodies data
    adata_abs = load_visium_with_abs_data()
    if adata_abs is not None and adata_abs.shape[0] > 0:
        adata_abs = basic_preprocessing(adata_abs)
        create_qc_plots(adata_abs, 'visium_with_abs')
        save_processed_data(adata_abs, 'visium_with_abs_processed.h5ad')
    else:
        print("⚠️ Skipping Visium+ABs data processing due to loading issues")
        adata_abs = None
    
    print("\n🎉 Data processing pipeline completed!")
    print("📁 Check the following directories:")
    print("   - processed_data/ : H5AD files ready for analysis")
    print("   - plots/ : Quality control plots")
    
    # Print summary
    if adata_spatial is not None:
        print(f"\n📊 Visium Spatial Summary:")
        print(f"   - Shape: {adata_spatial.shape}")
        print(f"   - Samples: {adata_spatial.obs['sample'].nunique() if 'sample' in adata_spatial.obs else 'N/A'}")
        print(f"   - Zones: {adata_spatial.obs['zonationGroup'].nunique() if 'zonationGroup' in adata_spatial.obs else 'N/A'}")
        
    if adata_abs is not None:
        print(f"\n📊 Visium + Antibodies Summary:")
        print(f"   - Shape: {adata_abs.shape}")
        print(f"   - Samples: {adata_abs.obs['orig.ident'].nunique() if 'orig.ident' in adata_abs.obs else 'N/A'}")
        print(f"   - Zones: {adata_abs.obs['zonationGroup'].nunique() if 'zonationGroup' in adata_abs.obs else 'N/A'}")
        print(f"   - Proteins: {adata_abs.obsm['protein'].shape[1] if 'protein' in adata_abs.obsm else 'N/A'}")

if __name__ == "__main__":
    main()
