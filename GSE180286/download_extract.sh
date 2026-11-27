#!/bin/bash

echo "========================================================================"
echo "DOWNLOADING GEO SUPPLEMENTARY DATA"
echo "========================================================================"

# Download the RAW tar file
echo ""
echo "Downloading GSE180286_RAW.tar..."
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180286/suppl/GSE180286_RAW.tar

if [ ! -f "GSE180286_RAW.tar" ]; then
    echo "Download failed!"
    exit 1
fi

echo "Download complete"

# Check file size
FILE_SIZE=$(du -h GSE180286_RAW.tar | cut -f1)
echo "File size: $FILE_SIZE"

# Create extraction directory
mkdir -p geo_raw_data
cd geo_raw_data

echo ""
echo "========================================================================"
echo "EXTRACTING TAR FILE"
echo "========================================================================"

# Extract
echo ""
echo "Extracting..."
tar -xvf ../GSE180286_RAW.tar

echo ""
echo "Extraction complete"

# List what we got
echo ""
echo "========================================================================"
echo "EXTRACTED FILES"
echo "========================================================================"
echo ""
ls -lh

# Count files
N_FILES=$(ls -1 | wc -l)
echo ""
echo "Total files: $N_FILES"

# Look for features files
echo ""
echo "Looking for features/genes files..."
ls -1 | grep -i "features\|genes" || echo "No features files found with standard names"

# Look for matrix files
echo ""
echo "Looking for matrix files..."
ls -1 | grep -i "matrix" || echo "No matrix files found"

# Look for barcodes
echo ""
echo "Looking for barcode files..."
ls -1 | grep -i "barcode" || echo "No barcode files found"

# Show sample of file types
echo ""
echo "File types present:"
ls -1 | sed 's/.*\.//' | sort | uniq -c

echo ""
echo "========================================================================"
echo "NEXT STEPS"
echo "========================================================================"
echo ""
echo "Files extracted to: geo_raw_data/"
echo ""
echo "Look for files named like:"
echo "  - *_features.tsv.gz (gene information)"
echo "  - *_matrix.mtx.gz (count matrix)"
echo "  - *_barcodes.tsv.gz (cell barcodes)"
echo ""
echo "Once you identify the features file, extract it:"
echo "  gunzip geo_raw_data/*features*.gz"
echo ""
echo "========================================================================"
