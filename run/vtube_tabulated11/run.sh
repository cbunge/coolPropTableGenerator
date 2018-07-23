#!/bin/bash
m4 system/circinlet.m4 > system/blockMeshDict
blockMesh
decomposePar
mpirun -np 48 tabularRhoCentralFoam -parallel
