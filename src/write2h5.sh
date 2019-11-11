#!/bin/bash

# Compiler converter
gfortran CBAmeshconverter.f -o converter.x -O3

# Ask user for filename
echo Enter filename in CBA format

read filename

# Run Converter
./converter.x <<< $filename


# Pack data into h5 file
python3 write2h5.py <<< $filename


# remove spare files
rm cellFace.txt
rm faceCell.txt
rm faceBC.txt
rm faceInfo.txt
rm faceNodes.txt
rm faceType.txt
rm cellType.txt
rm nodeVertex.txt
rm *.dat

echo Program Complete!
