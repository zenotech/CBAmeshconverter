# Load packages

import h5py
import numpy as np
import csv
import os

# Read in data

cellFace_data = np.loadtxt('cellFace.txt', dtype='i', skiprows=1)
cellType_data = np.loadtxt('cellType.txt', dtype='i', skiprows=1)
faceBC_data = np.loadtxt('faceBC.txt', dtype='i', skiprows=1)
faceCell_data = np.loadtxt('faceCell.txt', dtype='i', skiprows=1)
faceInfo_data = np.loadtxt('faceInfo.txt', dtype='i', skiprows=1)
faceNodes_data = np.loadtxt('faceNodes.txt', dtype='i', skiprows=1)
faceType_data = np.loadtxt('faceType.txt', dtype='i', skiprows=1)
nodeVertex_data = np.loadtxt('nodeVertex.txt', dtype='d', skiprows=1)

# Read in Attribute data
attribute_data = np.loadtxt('cellFace.txt', dtype='i', max_rows=1)

# Identify dataset size
l_cellFace=len(cellFace_data);
l_cellType=len(cellType_data);
l_faceBC=len(faceBC_data);
l_faceCell=len(faceCell_data);
l_faceInfo=len(faceInfo_data);
l_faceNodes=len(faceNodes_data);
l_faceType=len(faceType_data);
l_nodeVertex=len(nodeVertex_data);

# Create file type
name = input()
base=os.path.splitext(name)[0]
f = h5py.File(base + ".h5","w")

# Create groups

mesh = f.create_group("mesh")

# Assign attributes
mesh.attrs.create("numFaces", attribute_data[2], shape=(1,1))
mesh.attrs.create("numCells", attribute_data[1], shape=(1,1))

# Create datasets
cellFace = mesh.create_dataset("cellFace", data=cellFace_data)
cellType = mesh.create_dataset("cellType", data=cellType_data)
faceBC = mesh.create_dataset("faceBC", data=faceBC_data, shape=(l_faceBC,1))
faceCell = mesh.create_dataset("faceCell", data=faceCell_data)
faceInfo = mesh.create_dataset("faceInfo", data=faceInfo_data)
faceNodes = mesh.create_dataset("faceNodes", data=faceNodes_data, shape=(l_faceNodes,1))
faceType = mesh.create_dataset("faceType", data=faceType_data, shape=(l_faceType,1))
nodeVertex = mesh.create_dataset("nodeVertex", data=nodeVertex_data)
