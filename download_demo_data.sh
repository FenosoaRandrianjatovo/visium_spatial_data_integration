#!/bin/bash
# Download Demo Data for Visium Spatial Integration Pipeline
# ========================================================
# This script helps users download the required demo data files

set -e  # Exit on any error

echo "📥 VISIUM SPATIAL DATA DOWNLOAD SCRIPT"
echo "====================================="
echo ""

# Create required directory structure
echo "🏗️  Creating directory structure..."
mkdir -p Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium
mkdir -p Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium

echo "✅ Directory structure created"
echo ""

# Display required files
echo "📋 REQUIRED DATA FILES:"
echo "----------------------"
echo "You need to obtain these files from your data source:"
echo ""
echo "🔬 Standard Visium Spatial Data:"
echo "   Data/Visium_Spatial/annot_mouseStStVisium.csv"
echo "   Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz"
echo "   Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/barcodes.tsv.gz" 
echo "   Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/features.tsv.gz"
echo ""
echo "🧪 Visium + Antibodies Data:"
echo "   Data/Visium_spatial_with_ABs/annot_mouseStStWithABsVisium.csv"
echo "   Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_mouseStStVisiumWithABs.h5"
echo "   Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_ADT_mouseStStVisiumWithABs.csv"
echo ""

# Data source information
echo "🌐 DATA SOURCES:"
echo "---------------"
echo "These files can be obtained from:"
echo ""
echo "📊 10X Genomics Public Datasets:"
echo "   https://www.10xgenomics.com/resources/datasets"
echo ""
echo "📚 GEO (Gene Expression Omnibus):"
echo "   https://www.ncbi.nlm.nih.gov/geo/"
echo "   Search for: 'Visium spatial transcriptomics liver'"
echo ""
echo "🔗 Original Publication Data:"
echo "   Check the methods section of spatial transcriptomics papers"
echo "   Look for liver zonation or hepatic spatial studies"
echo ""

# File format requirements
echo "📝 FILE FORMAT REQUIREMENTS:"
echo "----------------------------"
echo ""
echo "🔢 10X Matrix Files (standard Cell Ranger output):"
echo "   - matrix.mtx.gz    : Sparse count matrix"
echo "   - barcodes.tsv.gz  : Spot/cell barcodes"
echo "   - features.tsv.gz  : Gene identifiers"
echo ""
echo "💽 H5 Files:"
echo "   - .h5 format with '/matrix' group containing:"
echo "     - data, indices, indptr (sparse matrix)"
echo "     - barcodes, features (identifiers)"
echo ""
echo "📋 Annotation Files (CSV format):"
echo "   Required columns:"
echo "   - 'spot' or 'barcode': Cell/spot identifiers"
echo "   - 'sample': Sample identifiers (JBO1, JBO2, etc.)"
echo "   - 'zonationGroup': Tissue zones (Periportal, Central, Mid, Portal)"
echo "   - Additional metadata columns are welcome"
echo ""

# Example download commands (commented out - users need real URLs)
echo "💡 EXAMPLE DOWNLOAD COMMANDS:"
echo "-----------------------------"
echo "Once you have the URLs, uncomment and modify these commands:"
echo ""
echo "# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz 'YOUR_MATRIX_URL'"
echo "# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/barcodes.tsv.gz 'YOUR_BARCODES_URL'"
echo "# wget -O Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/features.tsv.gz 'YOUR_FEATURES_URL'"
echo "# wget -O Data/Visium_Spatial/annot_mouseStStVisium.csv 'YOUR_ANNOTATION_URL'"
echo ""
echo "# curl -o Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_mouseStStVisiumWithABs.h5 'YOUR_H5_URL'"
echo "# curl -o Data/Visium_spatial_with_ABs/annot_mouseStStWithABsVisium.csv 'YOUR_AB_ANNOTATION_URL'"
echo ""

# Validation
echo "🔍 VALIDATION:"
echo "-------------"
echo "After downloading, validate your data structure:"
echo ""
echo "python test_pipeline.py"
echo ""

# Next steps
echo "🚀 NEXT STEPS:"
echo "-------------"
echo "1. Download the required files using the methods above"
echo "2. Ensure files are in the correct directory structure"
echo "3. Run: python test_pipeline.py (to validate)"
echo "4. Run: python run_pipeline.py (to execute pipeline)"
echo ""

echo "✅ Setup script completed!"
echo "📧 Need help? Check the README.md for detailed instructions"
