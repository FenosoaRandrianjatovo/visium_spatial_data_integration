# 📊 Data Structure Documentation

This document explains the required data structure and formats for the Visium Spatial Integration Pipeline.

## 📁 Required Directory Structure

```
Data/
├── Visium_Spatial/                                    # Standard Visium spatial transcriptomics
│   ├── annot_mouseStStVisium.csv                     # Spot annotations and metadata
│   └── rawData_mouseStStVisium/
│       ├── sampleComp_mouseStStVisium.txt            # Sample composition info (optional)
│       └── countTable_mouseStStVisium/               # 10X Cell Ranger output
│           ├── matrix.mtx.gz                         # Sparse count matrix
│           ├── barcodes.tsv.gz                       # Spot barcodes
│           └── features.tsv.gz                       # Gene features
│
└── Visium_spatial_with_ABs/                          # Visium + Antibody data
    ├── annot_mouseStStWithABsVisium.csv              # Spot annotations with protein info
    └── rawData_mouseStStWithABsVisium/
        ├── countTable_mouseStStVisiumWithABs.h5      # Combined RNA+protein data (H5 format)
        └── countTable_ADT_mouseStStVisiumWithABs.csv # Antibody-derived tag counts
```

## 📋 File Format Specifications

### **Annotation Files (CSV)**

#### `annot_mouseStStVisium.csv`
Required columns:
- `UMAP_1`, `UMAP_2`: Spatial coordinates
- `spot`: Spot barcode identifier  
- `sample`: Sample identifier (JBO1, JBO2, JBO3, JBO4, CAP02)
- `type`: Spot type ('tissue' or 'background')
- `cluster`: Cluster assignment (0, 1, 2, 3, 4)
- `zonation`: Continuous zonation score (0.0-1.0)
- `zonationGroup`: Discrete zone ('Periportal', 'Mid', 'Central', 'Portal')

Example:
```csv
UMAP_1,UMAP_2,spot,sample,type,cluster,zonation,zonationGroup
-1.70713826,1.90631304,AAACACCAATAACTGC-1_1,JBO1,tissue,4,0.50615025,Mid
2.96770152,1.06804667,AAACATTTCCCGGATT-1_1,JBO1,tissue,3,0.6382178,Mid
```

#### `annot_mouseStStWithABsVisium.csv`
Required columns:
- `spot`: Spot barcode identifier
- `orig.ident`: Sample identifier  
- `nCount_Spatial`, `nFeature_Spatial`: QC metrics
- `percent.mt`: Mitochondrial gene percentage
- `zonation`: Continuous zonation score
- `zonationGroup`: Discrete zone
- `CD4-adtSignal`, `CD8a-adtSignal`, etc.: Protein abundance values

### **Count Matrix Files**

#### **10X Format (`countTable_mouseStStVisium/`)**
Standard Cell Ranger output:
- `matrix.mtx.gz`: Market Matrix format, sparse count matrix
- `barcodes.tsv.gz`: One barcode per line (spot identifiers)
- `features.tsv.gz`: Gene information (may have 1-3 columns)

#### **H5 Format (`countTable_mouseStStVisiumWithABs.h5`)**
HDF5 file with `/matrix` group containing:
- `data`: Non-zero count values
- `indices`: Row indices for sparse matrix
- `indptr`: Column pointers for sparse matrix
- `shape`: Matrix dimensions
- `barcodes`: Spot identifiers
- `features`: Gene/feature information

#### **ADT Counts (`countTable_ADT_mouseStStVisiumWithABs.csv`)**
CSV file with:
- Rows: Spot barcodes (matching main dataset)
- Columns: Antibody markers (CD4, CD8a, CD117, etc.)
- Values: Protein abundance measurements

## 🔍 Data Characteristics

### **Spatial Dataset Stats**
- **Spots**: ~5,862 after QC filtering
- **Genes**: ~16,887 detected genes
- **Samples**: 5 biological samples
- **Zones**: 5 zonation groups (including background)
- **Technology**: 10X Visium spatial gene expression

### **Antibody Dataset Stats**  
- **Spots**: ~1,659 after QC filtering
- **Genes**: ~14,671 RNA features
- **Proteins**: ~80 antibody markers
- **Samples**: 1 sample with multi-modal data
- **Technology**: Visium + protein detection

## 🧬 Biological Context

### **Sample Information**
- **JBO1-JBO4**: Biological replicates of liver tissue
- **CAP02**: Additional liver sample
- **BH1**: Multi-modal sample with protein data

### **Liver Zonation**
- **Periportal**: Zone 1, near portal vein
- **Mid**: Intermediate metabolic zone  
- **Central**: Zone 3, near central vein
- **Portal**: Portal tract regions

### **Protein Markers**
Key antibodies for liver analysis:
- **Immune markers**: CD4, CD8a, CD19, CD3e
- **Metabolic markers**: CD71, CD326
- **Structural markers**: CD31, CD44

## ⚙️ Technical Requirements

### **File Size Expectations**
- `matrix.mtx.gz`: ~100-200 MB
- `annot_*.csv`: ~1-5 MB  
- `*.h5`: ~50-100 MB
- Total dataset: ~300-500 MB

### **Format Validation**
The pipeline automatically validates:
- File existence and readability
- Required columns in annotations
- Matrix format compatibility
- Barcode matching between files

### **Common Issues**
1. **Barcode mismatch**: Ensure annotation barcodes match matrix barcodes
2. **Missing columns**: Check required columns in annotation files
3. **File corruption**: Verify complete download of compressed files
4. **Format differences**: 10X output may vary between Cell Ranger versions

## 📥 Data Acquisition Guide

### **Public Data Sources**
1. **10X Genomics Datasets**: https://www.10xgenomics.com/resources/datasets
2. **GEO Database**: Search "Visium spatial liver zonation"
3. **SRA**: NCBI Sequence Read Archive
4. **Original Publications**: Check supplementary data

### **Data Processing Notes**
- Pipeline handles different 10X Cell Ranger versions
- Automatically detects feature file format (1-3 columns)
- Robust to minor annotation format variations
- Includes error handling for common data issues

## 🔧 Customization for Your Data

### **Adapting to Your Dataset**
1. **Match directory structure** exactly as shown
2. **Update column names** in annotation files if needed
3. **Modify sample identifiers** in scripts if using different names
4. **Adjust zonation labels** for your tissue type

### **Required Modifications**
If using different tissue types:
- Update `zonationGroup` values to match your biology
- Modify cluster interpretation in analysis scripts
- Adjust QC thresholds for your data characteristics

---

**📧 Questions about data format?** Check the main README.md or run `python test_pipeline.py` to validate your data structure.
