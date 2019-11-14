#!/bin/bash

# Request mesh file name

echo Enter .h5 mesh file name

read filename

# Remove existing dump files
rm cellFace.txt > /dev/null 2>&1
rm cellInfo.txt > /dev/null 2>&1
rm faceBC.txt > /dev/null 2>&1
rm faceCell.txt > /dev/null 2>&1
rm faceInfo.txt > /dev/null 2>&1
rm faceNodes.txt > /dev/null 2>&1
rm faceType.txt > /dev/null 2>&1
rm nodeVertex.txt > /dev/null 2>&1

# Dump out files to text
h5dump -d /mesh/cellFace $filename >> cellFace.txt
h5dump -d /mesh/cellType $filename >> cellInfo.txt
h5dump -d /mesh/faceBC $filename >> faceBC.txt
h5dump -d /mesh/faceCell $filename >> faceCell.txt
h5dump -d /mesh/faceInfo $filename >> faceInfo.txt
h5dump -d /mesh/faceNodes $filename >> faceNodes.txt
h5dump -d /mesh/faceType $filename >> faceType.txt
h5dump -d /mesh/nodeVertex $filename >> nodeVertex.txt

echo Program Complete!
