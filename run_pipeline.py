#!/usr/bin/env python3
"""
Master Pipeline Runner
=====================

This script runs the complete Visium spatial transcriptomics analysis pipeline.
It orchestrates data processing, spatial integration, and downstream analysis.
"""

import subprocess
import sys
import os
import time
from pathlib import Path

def run_command(command, description):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"🚀 {description}")
    print(f"{'='*60}")
    print(f"Running: {command}")
    
    start_time = time.time()
    
    try:
        # Use the virtual environment Python
        if command.startswith('python'):
            command = command.replace('python', '/lustre10/scratch/fenosoa/visium_spatial_integration/.venv/bin/python', 1)
        
        result = subprocess.run(command, shell=True, check=True, 
                              capture_output=True, text=True)
        
        elapsed_time = time.time() - start_time
        print(f"✅ {description} completed successfully in {elapsed_time:.1f}s")
        
        # Print last few lines of output for confirmation
        if result.stdout:
            output_lines = result.stdout.strip().split('\n')
            print("📋 Output summary:")
            for line in output_lines[-5:]:
                print(f"   {line}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print(f"❌ {description} failed after {elapsed_time:.1f}s")
        print(f"Error: {e}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        return False

def check_environment():
    """Check if the Python environment is properly set up."""
    print("🔍 Checking Python environment...")
    
    # Check if virtual environment exists
    venv_path = Path("/lustre10/scratch/fenosoa/visium_spatial_integration/.venv")
    if not venv_path.exists():
        print("❌ Virtual environment not found!")
        print("Let's run: python -m venv .venv")
        subprocess.run(
            "python -m venv .venv",
            shell=True, check=True, capture_output=True, text=True
        )
        print("✅ Virtual environment created")
        # return True
    
    # Check if required packages are installed
    try:
        result = subprocess.run(
            "/lustre10/scratch/fenosoa/visium_spatial_integration/.venv/bin/python -c 'import scanpy, pandas, anndata'",
            shell=True, check=True, capture_output=True, text=True
        )
        print("✅ Core packages available")
        return True
    except subprocess.CalledProcessError:
        print("❌ Required packages not installed!")
        print("Installing packages from requirements.txt...")
        return run_command(
            "/lustre10/scratch/fenosoa/visium_spatial_integration/.venv/bin/pip install -r requirements.txt",
            "Installing required packages"
        )

def check_data():
    """Check if the required data files exist."""
    print("🔍 Checking data files...")
    
    required_files = [
        "Data/Visium_Spatial/annot_mouseStStVisium.csv",
        "Data/Visium_Spatial/rawData_mouseStStVisium/countTable_mouseStStVisium/matrix.mtx.gz",
        "Data/Visium_spatial_with_ABs/annot_mouseStStWithABsVisium.csv",
        "Data/Visium_spatial_with_ABs/rawData_mouseStStWithABsVisium/countTable_mouseStStVisiumWithABs.h5"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print("❌ Missing data files:")
        for file_path in missing_files:
            print(f"   - {file_path}")
        return False
    else:
        print("✅ All required data files found")
        return True

def run_pipeline_step(script_name, description, required=True):
    """Run a pipeline step and return success status."""
    success = run_command(f"python {script_name}", description)
    
    if not success and required:
        print(f"\n❌ Pipeline failed at: {description}")
        print("Cannot continue to next steps.")
        return False
    elif not success:
        print(f"\n⚠️ Optional step failed: {description}")
        print("Continuing with remaining steps...")
    
    return success

def main():
    """Run the complete pipeline."""
    print("🚀 VISIUM SPATIAL TRANSCRIPTOMICS ANALYSIS PIPELINE")
    print("=" * 60)
    print("This pipeline will:")
    print("1. Process Visium spatial data")
    print("2. Integrate datasets using SpatialGlue")
    print("3. Perform downstream analysis")
    print("=" * 60)
    
    # Pre-flight checks
    if not check_environment():
        print("❌ Environment check failed. Please fix issues and try again.")
        return
    
    if not check_data():
        print("❌ Data check failed. Please ensure all data files are present.")
        return
    
    print("\n✅ All pre-flight checks passed!")
    
    # Pipeline steps
    pipeline_steps = [
        ("data_processing.py", "Data Processing and Quality Control", True),
        ("spatial_integration.py", "Spatial Data Integration", True),
        ("downstream_analysis.py", "Downstream Analysis", False)
    ]
    
    successful_steps = 0
    total_start_time = time.time()
    
    for script, description, required in pipeline_steps:
        if run_pipeline_step(script, description, required):
            successful_steps += 1
    
    # Pipeline summary
    total_elapsed = time.time() - total_start_time
    print(f"\n{'='*60}")
    print("🎉 PIPELINE SUMMARY")
    print(f"{'='*60}")
    print(f"✅ Completed steps: {successful_steps}/{len(pipeline_steps)}")
    print(f"⏱️ Total runtime: {total_elapsed/60:.1f} minutes")
    
    if successful_steps == len(pipeline_steps):
        print("\n🎉 All pipeline steps completed successfully!")
        print("\n📁 Output directories:")
        print("   - processed_data/    : Preprocessed individual datasets")
        print("   - integrated_data/   : Integrated spatial data")
        print("   - plots/             : Quality control and integration plots")
        print("   - results/           : Differential expression and spatial analysis")
        
        print("\n📊 Key files to check:")
        print("   - plots/visium_spatial_qc_metrics.png")
        print("   - plots/visium_with_abs_qc_metrics.png")
        print("   - plots/integration/spatialglue_integration_overview.png")
        print("   - results/summary/integration_summary.png")
        print("   - results/summary/analysis_summary.txt")
        
    else:
        print(f"\n⚠️ Pipeline completed with {len(pipeline_steps) - successful_steps} failures")
        print("Check the error messages above for troubleshooting.")
    
    print(f"\n{'='*60}")

if __name__ == "__main__":
    main()
