# 🧬 Visium Spatial Transcriptomics Integration Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Spatial Analysis](https://img.shields.io/badge/analysis-spatial-green.svg)](https://www.10xgenomics.com/products/spatial-gene-expression)

A complete end-to-end pipeline for processing and integrating **Visium spatial transcriptomics data**. This pipeline handles multi-modal data (RNA + Protein), performs quality control, spatial integration, and generates publication-ready visualizations.

## 🎯 Pipeline Overview

### **Why Spatial Data Integration?**

- **🔬 Multi-sample Analysis**: Integrate data across different samples/conditions
- **🧪 Multi-modal Integration**: Combine RNA and protein data from the same tissue  
- **📍 Spatial Context Preservation**: Maintain spatial relationships during integration
- **⚖️ Batch Effect Correction**: Remove technical artifacts while preserving biological variation
- **🗺️ Zonation Analysis**: Study tissue organization patterns (e.g., liver zonation)

### **Key Features**

✅ **Robust Data Loading**: Handles 10X matrix and H5 formats automatically  
✅ **Quality Control**: Comprehensive QC metrics and filtering  
✅ **Spatial Integration**: Multi-dataset harmonization with Harmony/SpatialGlue  
✅ **Visualization**: Publication-ready plots and spatial maps  
✅ **Reproducible**: Virtual environment ensures consistent results  

## 📊 Expected Results

This pipeline successfully processes:
- **Dataset 1**: Visium Spatial (5,862 spots × 16,887 genes, 5 samples)
- **Dataset 2**: Visium + Antibodies (1,659 spots × 14,671 genes, 80+ proteins)
- **Integration**: Combined 7,521 spots with preserved spatial organization

### **Generated Outputs**

1. **Quality Control Plots**: Gene/cell distributions, sample compositions
2. **Spatial Maps**: UMAP embeddings colored by zonation patterns  
3. **Integration Results**: Before/after comparison, batch correction assessment
4. **Biological Insights**: Liver zonation patterns (Periportal → Central → Portal)

## 🚀 Quick Start

### **1. Clone Repository**

```bash
git clone https://github.com/your-username/visium-spatial-integration.git
cd visium-spatial-integration
```

### **2. Prepare Input Data**

**⚠️ IMPORTANT**: The raw data files are too large for GitHub. You need to download them separately.

#### **Download the Demo Dataset**

This pipeline uses mouse liver Visium spatial transcriptomics data. Create the following directory structure:

```bash
mkdir -p Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium
mkdir -p Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium
```

#### **Required Files**

You need to obtain these files (see **Data Sources** section below):

```
Data/
├── Visium_Spatial/
│   ├── annot_mouseStStVisium.csv                    # Spot annotations
│   └── rawData_mouseStStVisium/
│       └── countTable_mouseStStVisium/
│           ├── barcodes.tsv.gz                      # 10X barcodes
│           ├── features.tsv.gz                      # 10X features  
│           └── matrix.mtx.gz                        # 10X count matrix
└── Visium_spatial_with_ABs/
    ├── annot_mouseStStWithABsVisium.csv            # Antibody data annotations
    └── rawData_mouseStStWithABsVisium/
        ├── countTable_ADT_mouseStStVisiumWithABs.csv    # ADT counts
        └── countTable_mouseStStVisiumWithABs.h5         # H5 format data
```

### **3. Run Pipeline**

```bash
# Run the complete analysis (auto-installs dependencies)
python run_pipeline.py
```

**That's it!** The pipeline will:
- Create virtual environment automatically
- Install all required packages  
- Process both datasets
- Perform spatial integration
- Generate all plots and results

### **4. View Results**

```bash
# Check generated outputs
ls plots/                    # Quality control and integration plots
ls processed_data/          # Clean H5AD files
ls integrated_data/         # Final integrated dataset
```

## 📁 Project Structure

```
visium-spatial-integration/
├── 📜 Core Pipeline Scripts
│   ├── data_processing.py           # Data loading & preprocessing
│   ├── spatial_integration.py       # Multi-dataset integration  
│   ├── downstream_analysis.py       # Advanced analysis
│   └── run_pipeline.py             # Master pipeline runner
│
├── 🔧 Utilities & Testing
│   ├── test_pipeline.py            # Validation tests
│   ├── showcase_results.py         # Results demonstration
│   └── pipeline_summary.py         # Summary statistics
│
├── 📋 Documentation  
│   ├── README.md                    # This file
│   ├── requirements.txt             # Python dependencies
│   └── PIPELINE_SUCCESS.md         # Detailed results
│
├── 📊 Input Data (you download)
│   └── Data/                        # Raw Visium data files
│
└── 📈 Generated Outputs (pipeline creates)
    ├── processed_data/              # Clean datasets
    ├── integrated_data/             # Integrated results
    └── plots/                       # All visualizations
```

## 📥 Data Sources & Preparation

### **Option 1: Use Demo Dataset (Recommended)**

The pipeline was designed and tested with mouse liver spatial transcriptomics data:

#### **Standard Visium Data**
- **Source**: 10X Genomics public datasets or GEO
- **Format**: 10X Cell Ranger output (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
- **Sample info**: Contains 5 biological samples (JBO1-JBO4, CAP02)
- **Annotations**: Includes liver zonation information (Periportal, Central, Mid, Portal)

#### **Visium + Antibodies Data**  
- **Source**: Multi-modal Visium dataset
- **Format**: H5 file + ADT CSV
- **Features**: ~80 protein markers alongside RNA
- **Spatial info**: Same tissue with protein quantification

### **Option 2: Use Your Own Data**

To use your own Visium data:

1. **Organize your files** to match the expected structure
2. **Update annotation column names** in the scripts if needed
3. **Modify sample identifiers** in the processing scripts

#### **File Format Requirements**

- **10X Matrix**: Standard Cell Ranger output
- **H5 Files**: 10X HDF5 format with `/matrix` group
- **Annotations**: CSV with columns:
  - `spot`: Barcode identifiers
  - `sample`: Sample/batch identifiers  
  - `zonationGroup`: Biological groupings (optional)
  - Additional metadata columns

### **Example Data Download Script**

Create this script to download demo data:

```bash
#!/bin/bash
# download_demo_data.sh

echo "📥 Downloading Visium demo data..."

# Create directories
mkdir -p Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium
mkdir -p Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium

# Download from your data source
# Replace these URLs with actual data sources
# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz "YOUR_MATRIX_URL"
# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/barcodes.tsv.gz "YOUR_BARCODES_URL"
# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/features.tsv.gz "YOUR_FEATURES_URL"

echo "✅ Data download complete!"
echo "🚀 Ready to run: python run_pipeline.py"
```

## 🔬 Scientific Applications

### **Demonstrated Applications**

- **🧬 Liver Zonation Analysis**: Spatial organization of hepatic functions
- **🔄 Multi-sample Integration**: Comparison across biological conditions
- **🧪 Multi-modal Analysis**: RNA + protein spatial co-localization  
- **📊 Quality Control**: Comprehensive data validation workflows

### **Potential Extensions**

- **Disease vs. Healthy**: Compare pathological tissue organization
- **Development Studies**: Track spatial patterns over time
- **Drug Effects**: Assess treatment impact on tissue architecture
- **Species Comparison**: Cross-species spatial organization

## ⚙️ Technical Details

### **Integration Methods**

- **Primary**: Harmony batch correction (robust, fast)
- **Alternative**: SpatialGlue (when available)
- **Fallback**: Scanpy basic integration

### **Key Parameters**

```python
# Preprocessing
min_genes = 200        # Minimum genes per spot
min_cells = 3         # Minimum spots per gene  
target_sum = 1e4      # Normalization target

# Integration  
n_neighbors = 15      # UMAP neighbors
n_pcs = 50           # Principal components
resolution = 0.5      # Clustering resolution
```

### **System Requirements**

- **Python**: 3.8+ (automatically managed via virtual environment)
- **Memory**: 8GB+ RAM recommended
- **Storage**: 2GB+ free space for outputs
- **OS**: macOS, Linux, Windows

## 🎯 Validation & Testing

### **Pipeline Validation**

```bash
# Test pipeline setup
python test_pipeline.py

# View pipeline capabilities  
python showcase_results.py

# Generate summary statistics
python pipeline_summary.py
```

### **Expected Runtime**

- **Data Loading**: ~2-3 minutes
- **Integration**: ~5-10 minutes  
- **Visualization**: ~2-3 minutes
- **Total**: ~10-15 minutes for complete analysis

## 🤝 Contributing

### **Adding New Features**

1. Fork the repository
2. Create feature branch: `git checkout -b feature-name`
3. Add your improvements
4. Test with `python test_pipeline.py`
5. Submit pull request

### **Reporting Issues**

Please include:
- Error messages and logs
- Data format details
- System information
- Steps to reproduce

## 📚 Dependencies

All automatically installed via `requirements.txt`:

```
scanpy>=1.9.0          # Single-cell analysis
pandas>=1.5.0          # Data manipulation  
anndata>=0.8.0         # Data structures
matplotlib>=3.5.0      # Plotting
seaborn>=0.11.0        # Statistical visualization
harmonypy>=0.0.6       # Batch correction
numpy<2.3              # Numerical computing (version locked)
```

## 📖 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{visium_spatial_integration,
  title={Visium Spatial Transcriptomics Integration Pipeline},
  author={Your Name},
  year={2025},
  url={https://github.com/your-username/visium-spatial-integration}
}
```

**Key Methods:**
- **Harmony**: Korsunsky et al. (2019) Nature Biotechnology
- **Scanpy**: Wolf et al. (2018) Genome Biology  
- **10X Visium**: 10X Genomics Spatial Gene Expression Solution

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 🏆 Success Stories

### **Validated Results**

✅ **5,862 spatial spots** successfully processed  
✅ **16,887 genes** analyzed across multiple samples  
✅ **Liver zonation patterns** clearly identified  
✅ **Multi-modal integration** of RNA + 80 proteins  
✅ **Publication-ready figures** generated automatically  

---

**🚀 Ready to explore spatial biology? Clone this repo and discover the spatial organization of your tissue!**

```bash
git clone https://github.com/your-username/visium-spatial-integration.git
cd visium-spatial-integration
python run_pipeline.py
```