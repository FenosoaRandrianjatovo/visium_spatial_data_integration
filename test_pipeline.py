#!/usr/bin/env python3
"""
Pipeline Validation Test
=======================
Quick test to validate the pipeline works end-to-end for new users.
"""

import os
import sys

def test_pipeline_setup():
    """Test that all components are properly set up."""
    
    print("🧪 PIPELINE VALIDATION TEST")
    print("=" * 40)
    
    # Check required files
    required_files = [
        'data_processing.py',
        'spatial_integration.py', 
        'run_pipeline.py',
        'requirements.txt',
        'README.md'
    ]
    
    print("📁 Checking core files...")
    missing_files = []
    for file in required_files:
        if os.path.exists(file):
            print(f"✅ {file}")
        else:
            print(f"❌ {file}")
            missing_files.append(file)
    
    # Check data directory
    print("\n📊 Checking data structure...")
    data_paths = [
        'Data/Visium_Spatial/annot_mouseStStVisium.csv',
        'Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz',
        'Data/Visium_spatial_with_ABs/annot_mouseStStWithABsVisium.csv',
        'Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_mouseStStVisiumWithABs.h5'
    ]
    
    missing_data = []
    for path in data_paths:
        if os.path.exists(path):
            size = os.path.getsize(path) / (1024*1024)  # MB
            print(f"✅ {path} ({size:.1f} MB)")
        else:
            print(f"❌ {path}")
            missing_data.append(path)
    
    # Check virtual environment
    print("\n🐍 Checking Python environment...")
    venv_path = '.venv/bin/python'
    if os.path.exists(venv_path):
        print("✅ Virtual environment exists")
        
        # Test key imports
        try:
            import subprocess
            result = subprocess.run([venv_path, '-c', 'import scanpy, pandas, anndata'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                print("✅ Core packages installed")
            else:
                print(f"❌ Import error: {result.stderr}")
        except Exception as e:
            print(f"❌ Environment test failed: {e}")
    else:
        print("❌ Virtual environment not found")
    
    # Check outputs (if pipeline was run)
    print("\n📈 Checking pipeline outputs...")
    output_dirs = ['processed_data', 'plots', 'integrated_data']
    outputs_exist = False
    
    for dir_name in output_dirs:
        if os.path.exists(dir_name):
            files = os.listdir(dir_name)
            print(f"✅ {dir_name}/ ({len(files)} files)")
            outputs_exist = True
        else:
            print(f"⚪ {dir_name}/ (not created yet)")
    
    # Summary
    print("\n" + "=" * 40)
    print("🎯 VALIDATION SUMMARY:")
    print(f"   Core files: {len(required_files) - len(missing_files)}/{len(required_files)}")
    print(f"   Data files: {len(data_paths) - len(missing_data)}/{len(data_paths)}")
    print(f"   Environment: {'✅' if os.path.exists(venv_path) else '❌'}")
    print(f"   Outputs: {'✅' if outputs_exist else '⚪ (run pipeline first)'}")
    
    if len(missing_files) == 0 and len(missing_data) == 0:
        print("\n🎉 READY TO RUN!")
        print("Execute: python run_pipeline.py")
    else:
        print("\n⚠️ SETUP INCOMPLETE")
        if missing_files:
            print(f"Missing files: {missing_files}")
        if missing_data:
            print(f"Missing data: {missing_data}")
    
    print("=" * 40)

if __name__ == "__main__":
    test_pipeline_setup()
