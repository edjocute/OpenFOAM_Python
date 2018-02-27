
# coding: utf-8

# In[2]:

# Usage: python EOFv.py
# Input: Reads in a text file with extension ".VSGR"
# Input: Format of input file as follows
# Input: STL start
# Input: <STL Filenames> <levels of refinement of each | optional>
# Input: STL end
# Input: Refinement start
# Input: <STL Filenames around which refinement needed> <mode | optional> <levels | optional>
# Input: Refinement end
# Input: If <levels of refinement of each> not specified, then defaults to (2 3)
# Input: If <mode> not specified, defaults to "inside"
# Input: If <levels> not specified, defaults to (1.00 3)
# Output: OpenFOAM files created in directory from which EOF.py was called
# Notes: blockMeshDict created such that there are 10 cells in Z direction. 
# Notes: Number of cells in X & Y computed to ensure aspect ratio as close to 1 as possible
# Notes: No warranty on results. Use at your own risk/discretion
# Notes: Code is free. Appreciate feedSouth/acknowledging when using it
# Created by: Venugopalan Raghavan

# Update 26 May 2015:
# (1) Multiple solids in one STL file issue corrected
# (2) Unnecessary braces in refinementRegions definition corrected
# (3) Missing "Bottom" patch defined for variables
# (4) Incorrect class type for U corrected

# Update 23 June 2015:
# (1) Incorrect reference for name corrected - referenced "parts" instead of "part"

## Setup basic stuff

##Modifications by Bharathi:
# 1. Copy the stl files to constant/triSurface at the end
# 2. The boundary names are changed from left, right, front, back to West, East, North, South
# 3. decomposeParDict: method is changed from hierarchy to scotch and cellDecomposition is included in dataFile option
# 4. meshQualityDict file and its contents are added and accordingly meshQualityControls section in snappyHexMeshDict is modified
# 5. The write format is changed from ascii to binary in controlDict, else it's causing issues while reconstructing mesh.
# 6. surfaceFeatureExtractDict and its contents are included
# 7. blockMeshDict - extents follow 17% blockage. 
#The base resolution that corresponds to level 0 in snappy mesh is set to 20 m. 
#The number of cells division in x, y is based on this resolution, defined as l0resol.
#In the vertical direction, grading resolution is introduced
# 8. xrange is changed to range
# 9. Snap is turned off
# 10. features are commented in casetellated mesh
# 11. Refinement region - min Z is changed to 0
# 12. Few changes in controlDict, fvSchemes and fvSolution

##Modifications by Eddie (Feb 2018)
# 1. Introduced graded meshing in x,y,z directions and removed snappy refinements
# 2. Added 0/alphat, additional transport properties
# 3. Modified P_rgh to type fixedFluxPressure


gHeight=0.0

header = list()
header.append("/*------------------------------------------------------------------------*\\\n")
header.append("|=========                 |                                               |\n")
header.append("|\\\\      /   F ield        | OpenFOAM: The Open Source CFD Toolbox         |\n")
header.append("| \\\\    /    O peration    | Version:  2.1.0                               |\n")
header.append("|  \\\\  /     A nd          | Web:      www.OpenFOAM.org                    |\n")
header.append("|   \\\\/      M anipulation |                                               |\n")
header.append("\*------------------------------------------------------------------------*/\n")
header.append("\n")
header.append("\n")
header.append("FoamFile\n")
header.append("{\n")
header.append("/* Exported using EOFv */\n")
header.append("\tversion\t2.1;\n")
header.append("\tformat\tascii;\n")


# In[3]:

blockMesh = ["\tclass\tdictionary;\n","\tlocation\tconstant;\n","\tobject\tblockMeshDict;\n","}\n\n","convertToMeters\t1;\n"]

surfaceFeature = ["\tclass\tdictionary;\n","\tlocation\tsystem;\n","\tobject\tsurfaceFeatureExtractDict;\n","}\n\n"]

snappyMesh = ["\tclass\tdictionary;\n","\tlocation\tsystem;\n","\tobject\tsnappyHexMeshDict;\n","}\n\n","\ncastellatedMesh\ttrue;","\nsnap\tfalse;","\naddLayers\tfalse;"]

snappyAdd = list()
snappyAdd.append("\nsnapControls")
snappyAdd.append("\n{")
snappyAdd.append("\n\tnSmoothPatch\t2;")
snappyAdd.append("\n\ttolerance\t4;")
snappyAdd.append("\n\tnSolveIter\t20;")
snappyAdd.append("\n\tnRelaxIter\t4;")
snappyAdd.append("\n\tnFeatureSnapIter\t10;")
snappyAdd.append("\n}")
snappyAdd.append("\naddLayersControls")
snappyAdd.append("\n{")
snappyAdd.append("\n\trelativeSizes\ttrue;")
snappyAdd.append("\n\texpansionRatio\t1.3;")
snappyAdd.append("\n\tfinalLayerThickness\t0.4;")
snappyAdd.append("\n\tminThickness\t0.3;")
snappyAdd.append("\n\tnGrow\t0;")
snappyAdd.append("\n\tfeatureAngle\t45;")
snappyAdd.append("\n\tnRelaxIter\t4;")
snappyAdd.append("\n\tnSmoothSurfaceNormals\t1;")
snappyAdd.append("\n\tnSmoothNormals\t1;")
snappyAdd.append("\n\tnSmoothThickness\t10;")
snappyAdd.append("\n\tmaxFaceThicknessRatio\t0.4;")
snappyAdd.append("\n\tmaxThicknesstoMedialRatio\t4;")
snappyAdd.append("\n\tminMedianAxisAngle\t130;")
snappyAdd.append("\n\tnBufferCellsNoExtrude\t0;")
snappyAdd.append("\n\tnLayerIter\t30;")
snappyAdd.append("\n\tlayers")
snappyAdd.append("\n\t{")

snappyQual = list()
snappyQual.append("\n\nmeshQualityControls")
snappyQual.append("\n{")
snappyQual.append('\n\t#include\t"meshQualityDict"')
snappyQual.append("\n\tnSmoothScale\t4;")
snappyQual.append("\n\terrorReduction\t0.75;")
snappyQual.append("\n}")


# In[4]:

import numpy as np
import math

def nofcells(a,b,s):
    #ratio of last cell to first cell
    r=np.divide(b,a,dtype=np.float)
    k=np.divide((a-s),(b-s),dtype=np.float)
    n=1+np.divide(math.log(r),math.log(k),dtype=np.float)
    n1=np.rint(n)
    return n1

def grading_output(v1,v2,v3,v4,a,b):
    #Calculate the extents over which grading has to be performed. 
    #It is assumed here that there will be three sets of grading in the following order
    #Coarse to fine; Uniform; Fine to Coarse
    s1=np.fabs(v1-v2)
    s2=np.fabs(v3-v2)
    s3=np.fabs(v4-v3)
    #print(s1,s2,s3)

    #No. of cells and expansion ratios
    ncells=np.zeros(3,dtype='int')
    ncells[0]=nofcells(a,b,s1)
    ncells[1]=np.rint(np.divide(s2,b))
    ncells[2]=nofcells(b,a,s3)
    er_1=np.divide(b,a,dtype=np.float)
    er_2=np.divide(b,b)
    er_3=np.divide(a,b,dtype=np.float)

    cf=[s1,ncells[0],er_1]
    uni=[s2,ncells[1],er_2]
    fc=[s3,ncells[2],er_3]
    dirn=np.array([cf,uni,fc])

    return dirn[ncells>0],ncells.sum()


# In[5]:

import shutil
import os as os
cwd = os.getcwd().replace("\\","/")

os.mkdir("0")
os.mkdir("constant")
os.mkdir("constant/polyMesh")
os.mkdir("constant/triSurface")
os.mkdir("system")

iF = cwd+"/Input.VSGR"
iFR = open(iF,"r")
currline = iFR.readline()

STLFileList = list()
RefinementList = list()


# In[6]:

## Find out the names of the STL files and refinement regions and the levels of refinement

while currline!='':
	if "STL start" in currline:
		currline = iFR.readline()
		while "STL end" not in currline:
			name = currline.strip().replace('\t',' ')
			#dummy = dummy.split(None,1)
			print(name)
			#name = dummy[0]
			refSurf = "1"
			#if len(dummy) == 2:
				#assert dummy[1].startswith('(') and dummy[1].endswith(')'),'"'+dummy[1]+'" is not a valid refinement specification for '+name
			#	refSurf = dummy[1]
			#else:	
			#	refSurf = "1"
			STLFileList.append(name + "\t" + refSurf)
			currline = iFR.readline()
	elif "Refinement start" in currline:
		currline = iFR.readline()
		while "Refinement end" not in currline:
			dummy = currline.strip().split()
			name = dummy[0]
			if len(dummy) == 6:
				mode = dummy[1]
				refVol = ' '.join(dummy[2:])
			elif len(dummy) == 3:
				mode = "inside"
				refVol = ' '.join(dummy[1:])
			else:
				mode = "inside"
				refVol = "((1.00 0))"
			RefinementList.append(name + "\t" + mode + "\t" + refVol)
			currline = iFR.readline()
	currline = iFR.readline()

iFR.close()


# In[7]:

## Declare variables for the global maximum and minimum coordinates

gXMin = gYMin = gZMin = 1e10
gXMax = gYMax = gZMax = -1e10

refXMin = refYMin = refZMin = 1e10
refXMax = refYMax = refZMax = -1e10

## Find out the names of the solids in each STL file and the extents and write to MasterList

oF = cwd + "/MasterSTLList"
oFW = open(oF,"w")

for ii in range(0,len(STLFileList)):
	solidNames = list()
	localMinBounds = list()
	localMaxBounds = list()
	parts = STLFileList[ii].split("\t")
	name = parts[0]
	refSurf = parts[1]
	iF = cwd+"/"+name
	iFR = open(iF,"r")
	currline = iFR.readline()
	while currline!='':
		if ("solid" in currline) and ("endsolid" not in currline):
			xMin = yMin = zMin = 1e10
			xMax = yMax = zMax = -1e10
			solidNames.append(currline.strip().split()[1])
			currline = iFR.readline()
			while "endsolid" not in currline:		
				if "vertex" in currline:
					dummy = currline.strip().split()
					x = float(dummy[1])
					y = float(dummy[2])
					z = float(dummy[3])
					xMin = min(xMin,x)
					xMax = max(xMax,x)
					yMin = min(yMin,y)
					yMax = max(yMax,y)
					zMin = min(zMin,z)
					zMax = max(zMax,z)
				currline = iFR.readline()
			localMinBounds.append([xMin,yMin,zMin])
			localMaxBounds.append([xMax,yMax,zMax])
			gXMin = min(gXMin,xMin)
			gXMax = max(gXMax,xMax)
			gYMin = min(gYMin,yMin)
			gYMax = max(gYMax,yMax)
			gZMin = min(gZMin,zMin)
			gZMax = max(gZMax,zMax)
		currline = iFR.readline()
	iFR.close()
	for qq in range(0,len(solidNames)):
		oFW.write(name+"\t"+"\t"+solidNames[qq]+"\t"+str(localMinBounds)+"\t"+str(localMaxBounds)+"\n")

oFW.close()
iFR.close()

gmin_xy=([gXMin,gYMin])
gmax_xy=([gXMax,gYMax])
min_xyz=np.rint([gXMin,gYMin,0])
max_xyz=np.rint([gXMax,gYMax,gZMax])
diff_xyz=max_xyz-min_xyz
print "Min and Max vertices of STL regions are:"
print min_xyz
print max_xyz
print "The STL region extents are: " +str(diff_xyz)


## Set up the numbers for blockMeshDict

#computational domain size extents that satisfies 17% blockage (for both windward frontal and leeward sides)
#and domain height is 10 times the max value of z
#For 17%blockage, c=2.5
c=2.5
gXMin=min_xyz[0]-diff_xyz[0]*c
gXMax=max_xyz[0]+diff_xyz[0]*c
gYMin=min_xyz[1]-diff_xyz[1]*c
gYMax=max_xyz[1]+diff_xyz[1]*c
gZMin=min_xyz[2]
gZMax=8*diff_xyz[2]

#Domain lengths
Lx=gXMax-gXMin
Ly=gYMax-gYMin
Lz=gZMax
avgX = (gXMin + gXMax)/2
avgY = (gYMin + gYMax)/2
avgZ = (gZMin + gZMax)/2
print "The computational domain extents are:"
print(Lx,Ly,Lz)

#The min and max of computational domain
print "Min and Max vertices of comp. domain"
print (gXMin,gYMin,gZMin)
print (gXMax,gYMax,gZMax)

##### Blockmesh grading calculation

#first and last cell size and distance between last and first cell in horizontal plane
fine=2 #<------ USER INPUT HERE
coarse=fine*64
#Expand boundary by distance d in the horizontal plane 
d=0
xdirn,nx=grading_output(gXMin,min_xyz[0]-d,max_xyz[0]+d,gXMax,coarse,fine)
ydirn,ny=grading_output(gYMin,min_xyz[1]-d,max_xyz[1]+d,gYMax,coarse,fine)

#In z-direction, use uniform for first fz*gZMax and then graded after
finez=1
coarsez=64
fz = 0.02
zdirn,nz=grading_output(0,0,np.rint(fz*gZMax),gZMax,coarsez,finez)
#print "The total no. of cells:" +str(nz)
#print "Expansion ratio:" +str(er)

oF = cwd + "/constant/polyMesh/blockMeshDict"
oFW = open(oF,"w")

for ii in range(0,len(header)):
	oFW.write(header[ii])

for ii in range(0,len(blockMesh)):
	oFW.write(blockMesh[ii])

oFW.write("\n\nvertices\n(")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMin) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMin) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMax) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMax) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMin) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMin) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMax) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMax) + "\t" + str(gZMax) + ")")
oFW.write("\n);")
oFW.write("\nblocks\n(\n\thex\t(0 1 2 3 4 5 6 7)");
#oFW.write("("+str(dX)+" "+str(dY)+" "+str(dZ)+")"+"\tsimpleGrading\t(1 1 1)\n);") #uniform mesh
oFW.write("\t("+str(int(nx))+" "+str(int(ny))+" "+str(int(nz))+")")
oFW.write("\n\n\tsimpleGrading\n\t(")
for dirn in [xdirn,ydirn,zdirn]:
   oFW.write("\n\t\t(")
   for i in xrange(len(dirn)):
      oFW.write("\n\t\t\t({0[0]}\t{0[1]}\t{0[2]})".format(dirn[i]))
   oFW.write("\n\t\t)")
oFW.write("\n\t)\n);")
oFW.write("\nedges\n(\n);")
oFW.write("\nboundary\n(")
oFW.write("\n\tWest\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(0 4 7 3)\n\t\t);\n\t}")
oFW.write("\n\tEast\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(1 2 6 5)\n\t\t);\n\t}")
oFW.write("\n\tBottom\n\t{\n\t\ttype\twall;\n\t\tfaces\n\t\t(\n\t\t\t(0 3 2 1)\n\t\t);\n\t}")
oFW.write("\n\tTop\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(4 5 6 7)\n\t\t);\n\t}")
oFW.write("\n\tSouth\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(0 1 5 4)\n\t\t);\n\t}")
oFW.write("\n\tNorth\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(3 7 6 2)\n\t\t);\n\t}")
oFW.write("\n);")
oFW.write("\nmergePatchPairs\n(\n);")
oFW.close()


# In[8]:

## Setup and write snappyHexMeshDict

oF = cwd + "/system/snappyHexMeshDict" 
oFW = open(oF,"w")

for ii in range(0,len(header)):
	oFW.write(header[ii])

for ii in range(0,len(snappyMesh)):
	oFW.write(snappyMesh[ii])

oFW.write("\n\ngeometry\n{\n")

iF = cwd + "/MasterSTLList"
iFR = open(iF,"r")
currline = iFR.readline()

namesSoFar = list()

ii = 0

wSTL = list()
wSolid = list()
wRefSurf = list()
wMinBd = list()
wMaxBd = list()

while currline!='':
	part = currline.strip().split("\t")
	name = (part[0].split(".stl"))[0]
	refSurf = part[1]
	solid = part[2]
	minBd = part[3]
	maxBd = part[4]
	wSTL.append(name)
	wSolid.append(solid)
	wRefSurf.append(refSurf)
	wMinBd.append(minBd)
	wMaxBd.append(maxBd)
	currline = iFR.readline()

c1 = 0

for ii in range(0,len(wSTL)):
	name = wSTL[ii]
	refSurf = wRefSurf[ii]
	solid = wSolid[ii]
	minBd = wMinBd[ii]
	maxBd = wMaxBd[ii]
	if name not in namesSoFar:
		if ii!=0:
			#oFW.write("\n\t\t}\n\t}\n");
			oFW.write("\n\t}\n");
		namesSoFar.append(name)
		oFW.write("\n\t"+ name + ".stl\n\t{")
		oFW.write("\n\t\ttype\ttriSurfaceMesh;")
		oFW.write("\n\t\tname\t"+name+"_"+solid+";")
	else:
		oFW.write("\n\t\t\t"+solid)
		oFW.write("\n\t\t\t{")
		#oFW.write("\n\t\t\t\tname\t"+name+"_"+solid+";"+"\n\t\t\t}")
		oFW.write("\n\t\t\t\tname\t"+name+"_"+solid+";")
	if ii==len(wSTL)-1:
		#oFW.write("\n\t\t}\n\t}\n");
		oFW.write("\n\t}\n");

print(wSTL)

offset = 10

for ii in range(0,len(RefinementList)):
	parts = RefinementList[ii]
	name = (parts.split("\t")[0]).split(".stl")[0]
	pos = wSTL.index(name)
	minBd = wMinBd[pos]
	maxBd = wMaxBd[pos]
	verticesMin = minBd.split()
	verticesMax = maxBd.split()
	minX = (verticesMin[0].split(",")[0]).split("[[")[1] 
	minX = str(float(minX) - offset)    
	minY = (verticesMin[1].split(","))[0] 
	minY = str(float(minY) - offset)
	minZ = str(0) #minZ = (verticesMin[2].split("]"))[0] - minimumZ is set to 0
	maxX = (verticesMax[0].split(",")[0]).split("[[")[1]
	maxX = str(float(maxX) + offset)
	maxY = (verticesMax[1].split(","))[0]
	maxY = str(float(maxY) + offset)
	maxZ = (verticesMax[2].split("]"))[0]
	maxZ = str(float(maxZ) + offset)
	oFW.write("\n\t"+ name + "_Ref\n\t{")
	oFW.write("\n\t\ttype\tsearchableBox;")
	oFW.write("\n\t\tmin\t("+minX+" "+minY+" "+minZ+");")
	oFW.write("\n\t\tmax\t("+maxX+" "+maxY+" "+maxZ+");")
	oFW.write("\n\t}")

oFW.write("\n}")

oFW.write("\ncastellatedMeshControls")
oFW.write("\n{")
oFW.write("\n\tlocationInMesh\t("+str(avgX)+" "+str(avgY)+" "+str(avgZ)+");")
oFW.write("\n\tmaxLocalCells\t6000000;")
oFW.write("\n\tmaxGlobalCells\t60000000;")
oFW.write("\n\tminRefinementCells\t50;")
oFW.write("\n\tnCellsBetweenLevels\t3;")
oFW.write("\n\tresolveFeatureAngle\t60;")
oFW.write("\n\tallowFreeStandingZoneFaces\tfalse;")
oFW.write("\n\tfeatures")
oFW.write("\n\t(")
oFW.write("\n\t/*")
for ii in range(0,len(wSTL)):
	oFW.write("\n\t\t{")
	oFW.write('\n\t\t\tfile\t"'+wSTL[ii]+'.eMesh";')
	oFW.write('\n\t\t\tlevel\t2;')
	oFW.write("\n\t\t}")

oFW.write("\n\t\t*/")
oFW.write("\n\t);")

oFW.write("\n\trefinementSurfaces")
oFW.write("\n\t{")
namesSoFar = list()
for ii in range(0,len(wSTL)):
	name = wSTL[ii]
	#refSurf = wRefSurf[ii]
	refSurf="(0 0)"
	solid = wSolid[ii]
	if name not in namesSoFar:
		if ii!=0:
			oFW.write("\n\t\t}\n")
		#oFW.write("\n\t\t"+name)
                oFW.write("\n\t\t"+name+"_"+solid) #BB
		oFW.write("\n\t\t{")
		#oFW.write("\n\t\t\tlevel\t"+refSurf+";")
		#oFW.write("\n\t\t\tregions")
		#oFW.write("\n\t\t\t{")
		#oFW.write("\n\t\t\t\t"+name+"_"+solid)
		#oFW.write("\n\t\t\t\t{")
		#oFW.write("\n\t\t\t\t\tlevel\t"+refSurf+";")
		oFW.write("\n\t\t\tlevel\t"+refSurf+";")#BB
		oFW.write("\n\t\t\tpatchInfo")
		oFW.write("\n\t\t\t{")
		oFW.write("\n\t\t\t\ttype\twall;")
		oFW.write("\n\t\t\t}")
		#oFW.write("\n\t\t\t\t}")
	else:
		#oFW.write("\n\t\t\t\t"+name+"_"+solid)
		#oFW.write("\n\t\t\t\t{")
		#oFW.write("\n\t\t\t\t\tlevel\t"+refSurf+";")
		oFW.write("\n\t\t\tlevel\t"+refSurf+";") #BB
		oFW.write("\n\t\t\tpatchInfo")
		oFW.write("\n\t\t\t{")
		oFW.write("\n\t\t\t\ttype\twall;")
		oFW.write("\n\t\t\t}")
		#oFW.write("\n\t\t\t\t}")
	if ii==len(wSTL)-1:
		#oFW.write("\n\t\t\t}\n\t\t}\n")
		oFW.write("\n\t\t}\n")
oFW.write("\n\t}")

oFW.write("\n\trefinementRegions")
oFW.write("\n\t{")
#namesSoFar = list()
#for ii in range(0,len(RefinementList)):
#	parts = RefinementList[ii]
#	name = (parts.split("\t")[0]).split(".stl")[0]
#	mode = parts.split("\t")[1]
#	refVol = parts.split("\t")[2]
#	oFW.write("\n\t\t"+ name + "_Ref\n\t\t{")
#	oFW.write("\n\t\t\tmode\t"+mode+";")
#	oFW.write("\n\t\t\tlevels\t"+refVol+";")
#	oFW.write("\n\t\t}")

oFW.write("\n\t}")
oFW.write("\n}")

for ii in range(0,len(snappyAdd)):
	oFW.write(snappyAdd[ii])

for ii in range(0,len(wSTL)):
	oFW.write("\n\t\t"+wSTL[ii])
	oFW.write("\n\t\t{")
	oFW.write("\n\t\t\tnSurfaceLayers\t2;")
	oFW.write("\n\t\t}")

oFW.write("\n\t}")
oFW.write("\n}")

for ii in range(0,len(snappyQual)):
	oFW.write(snappyQual[ii])

oFW.write("\ndebug\t0;")
oFW.write("\nmergeTolerance\t1E-06;")
oFW.close()


# In[9]:

## Setup and write surfaceFeatureExtractDict
oF = cwd + "/system/surfaceFeatureExtractDict" 
oFW = open(oF,"w")

for ii in range(0,len(header)):
	oFW.write(header[ii])

for ii in range(0,len(surfaceFeature)):
	oFW.write(surfaceFeature[ii])

iF = cwd + "/MasterSTLList"
iFR = open(iF,"r")
currline = iFR.readline()

namesSoFar = list()

ii = 0

wSTL = list()
while currline!='':
	part = currline.strip().split("\t")
	name = (part[0].split(".stl"))[0]
	wSTL.append(name)
	currline = iFR.readline()

c1 = 0

for ii in range(0,len(wSTL)):
    name = wSTL[ii]
    if name not in namesSoFar:
        if ii!=0:
            oFW.write("\n\t}\n");
        namesSoFar.append(name)
        oFW.write("\n"+name+".stl\n\t{")
        oFW.write("\n\t\textractionMethod\textractFromSurface;")
        oFW.write("\n\t\textractFromSurfaceCoeffs")
        oFW.write("\n\t\t{")
        oFW.write("\n\t\t\tincludedAngle\t150;")
        oFW.write("\n\t\t}")
        oFW.write("\n\t\tsubsetFeatures")
        oFW.write("\n\t\t{")
        oFW.write("\n\t\t\tnonManifoldEdges\tno;")
        oFW.write("\n\t\t\topenEdges\tyes;")
        oFW.write("\n\t\t}")
    if ii==len(wSTL)-1:
        oFW.write("\n\t}\n")
                
oFW.close()


# In[12]:

## Create controlDict, decomposeParDict, fvSchemes, fvSolution, meshQualityDict and cuttingPlane

files = ["controlDict","decomposeParDict","fvSchemes","fvSolution","meshQualityDict","cuttingPlane"]

for f in files:
	oF = cwd + "/system/" + f
	oFW = open(oF,"w")
	
	for ii in range(0,len(header)):
		oFW.write(header[ii])
	
	oFW.write("\tclass\tdictionary;")
	oFW.write("\n\tlocation\tsystem;")
	oFW.write("\n\tobject\t"+f+";")
	oFW.write("\n}\n")

	if f == "controlDict":
		oFW.write("\napplication\tsimpleFoam;")
		oFW.write("\nstartFrom\tlatestTime;")
		oFW.write("\nstartTime\t0;")
		oFW.write("\nstopAt\tendTime;")
		oFW.write("\nendTime\t5000;")
		oFW.write("\ndeltaT\t1;")
		oFW.write("\nwriteControl\ttimeStep;")
		oFW.write("\nwriteInterval\t500;")
		oFW.write("\npurgeWrite\t1;")
		oFW.write("\nwriteFormat\tbinary;")
		oFW.write("\nwritePrecision\t6;")
		oFW.write("\nwriteCompression\tcompressed;")
		oFW.write("\ntimeFormat\tgeneral;")
		oFW.write("\ntimePrecision\t6;")
		oFW.write("\nrunTimeModifiable\ttrue;")
		oFW.write("\nfunctions")
		oFW.write("\n{")
		oFW.write("\n\t//#include \"probes\"")
		oFW.write("\n\t#include \"cuttingPlane\"")
		oFW.write("\n}")
	elif f == "decomposeParDict":
		oFW.write("\nnumberOfSubdomains\t1;")
#		oFW.write("\n\nmethod\thierarchical;\n")
		oFW.write("\n\nmethod\tscotch;\n")
		oFW.write("\nsimpleCoeffs")
		oFW.write("\n{")
		oFW.write("\n\tn\t(1 1 1);")
		oFW.write("\n\tdelta\t0.001;")
		oFW.write("\n}\n")
		oFW.write("\nhierarchicalCoeffs")
		oFW.write("\n{")
		oFW.write("\n\tn\t(1 1 1);")
		oFW.write("\n\tdelta\t0.001;")
		oFW.write("\n\torder\txyz;")
		oFW.write("\n}\n")
		oFW.write("\nmanualCoeffs")
		oFW.write("\n{")
		oFW.write('\n\tdataFile\t"cellDecomposition";')
		oFW.write("\n}")
	elif f == "fvSchemes":
		oFW.write("\n\nddtSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tsteadyState;")
		oFW.write("\n}")

		oFW.write("\n\ngradSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tGauss linear;")
		oFW.write("\n}")

		oFW.write("\n\ndivSchemes")
		oFW.write("\n{")	
		oFW.write("\n\tdefault\tnone;")
		oFW.write("\n\tdiv(phi,U)\tbounded Gauss upwind;")
		oFW.write("\n\tdiv(phi,T)\tbounded Gauss upwind;")
		oFW.write("\n\tdiv(phi,k)\tbounded Gauss upwind;")
		oFW.write("\n\tdiv(phi,epsilon)\tbounded Gauss upwind;")
		oFW.write("\n\tdiv(phi,omega)\tbounded Gauss upwind;")
		oFW.write("\n\tdiv((nuEff*dev(T(grad(U)))))\tGauss linear;")
		oFW.write("\n}")

		oFW.write("\n\nlaplacianSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tGauss linear limited 1;") #BB -commented below lines
		#oFW.write("\n\tdefault\tnone;")
		#oFW.write("\n\tlaplacian(nuEff,U)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian(kappaEff,T)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian(DkEff,k)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian(DepsilonEff,epsilon)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian(DomegaEff,omega)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian((1|A(U)),p)\tGauss linear limited 0.5;")
		#oFW.write("\n\tlaplacian((1|A(U)),p_rgh)\tGauss linear limited 0.5;")
		oFW.write("\n}")

		oFW.write("\n\ninterpolationSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tlinear;")
		oFW.write("\n}")

		oFW.write("\n\nsnGradSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tlimited 1;") #BB -commented below line
                #oFW.write("\n\tdefault\tlimited 0.5;")
		oFW.write("\n}")

		oFW.write("\n\nfluxRequired")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tno;")
		oFW.write("\n\tp	;")
		oFW.write("\n\tp_rgh	;")
		oFW.write("\n}")
	
	elif f == "fvSolution":
		oFW.write("\n\nsolvers")
		oFW.write("\n{")
		oFW.write('\n\t"p|p_rgh"')
		oFW.write("\n\t{")
		oFW.write("\n\t\tsolver\tGAMG;")
		oFW.write("\n\t\tsmoother\tGaussSeidel;")
		oFW.write("\n\t\ttolerance\t1e-05;")
		oFW.write("\n\t\trelTol\t0.01;")
		oFW.write("\n\t\tnPreSweeps\t0;")
		oFW.write("\n\t\tnPostSweeps\t2;")
		oFW.write("\n\t\tcacheAgglomeration\ton;")
		oFW.write("\n\t\tnCellsInCoarsestLevel\t10;")
		oFW.write("\n\t\tagglomerator\tfaceAreaPair;")
		oFW.write("\n\t\tmergeLevels\t1;")
		oFW.write("\n\t}")

		oFW.write('\n"k|omega|epsilon|U|T"')
		oFW.write("\n\t{")
		oFW.write("\n\t\tsolver\tsmoothSolver;")
		oFW.write("\n\t\tsmoother\tGaussSeidel;")
		oFW.write("\n\t\ttolerance\t1e-05;")
		oFW.write("\n\t\trelTol\t0.01;")
		oFW.write("\n\t}")
		oFW.write("\n}")

		oFW.write("\nSIMPLE")
		oFW.write("\n{")
		oFW.write("\n\tnNonOrthogonalCorrectors\t0;")
		oFW.write("\n\tpRefCell\t0;")
		oFW.write("\n\tpRefValue\t0;")
		oFW.write("\n\tresidualControl")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"k|omega|epsilon|U|T"\t1e-5;')
		oFW.write('\n\t\t"p|p_rgh"\t1e-5;')
		oFW.write("\n\t}")
		oFW.write("\n}")

		oFW.write("\nrelaxationFactors")
		oFW.write("\n{")
		oFW.write("\n\tequations")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"k|omega|epsilon|U|T"\t0.3;')
		oFW.write("\n\t}")
	
		oFW.write("\n\tfields")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"p|p_rgh"\t0.7;')
		oFW.write("\n\t}")
		
		oFW.write("\n}")

	elif f == "meshQualityDict":
		oFW.write("\n\n//\tInclude defaults parameters from master dictionary")
		oFW.write('\n#include "$WM_PROJECT_DIR/etc/caseDicts/meshQualityDict"')
		oFW.write("\n\n//\tminFaceWeight\t(0 -> 0.5)")
		oFW.write("\n\tminFaceWeight\t0.02;")
        
    	elif f == "cuttingPlane":
        	oFW.write("\n\ncuttingPlane")
        	oFW.write("\n\t{")
        	oFW.write("\n\t\ttype\t\t\tsurfaces;")
        	oFW.write('\n\t\tfunctionObjectLibs\t("libsampling.so");')
        	oFW.write("\n\t\toutputControl\toutputTime;")
        	oFW.write("\n\n\t\tsurfaceFormat\tvtk;")
        	oFW.write("\n\t\tfields \t\t( p U k epsilon T);")

        	oFW.write("\n\t\tinterpolationScheme\t\tcellPoint;")

        	oFW.write("\n\t\tsurfaces")
        	oFW.write("\n\t\t(")
        	oFW.write("\n\t\t\tzNormal")
        	oFW.write("\n\t\t\t{")
        	oFW.write("\n\t\t\t\ttype\t\t\tcuttingPlane;")
        	oFW.write("\n\t\t\t\tplaneType\t\t\tpointAndNormal;")
        	oFW.write("\n\t\t\t\tpointAndNormalDict")
        	oFW.write("\n\t\t\t\t{")
        	#oFW.write("\n\t\t\t\t\tbasePoint\t\t\t("+str(avgX)+" "+str(avgY)+" "+str(avgZ)+");")
                oFW.write("\n\t\t\t\t\tbasePoint\t\t\t("+str(avgX)+" "+str(avgY)+" "+str(20)+");") #BB
        	oFW.write("\n\t\t\t\t\tnormalVector\t\t\t(0 0 1);")
        	oFW.write("\n\t\t\t\t\t}")
        	oFW.write("\n\t\t\t\t\tinterpolate\t\t\ttrue;")
        	oFW.write("\n\t\t\t\t}")
        	oFW.write("\n\t\t\t);")
        	oFW.write("\n\t\t}")
    
	
	oFW.close()


# In[ ]:

## Create a changeDictionaryDict framework

oF = cwd + "/system/changeDictionaryDict"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tdictionary;")
oFW.write("\n\tlocation\tsystem;")
oFW.write("\n\tobject\tchangeDictionaryDict;")
oFW.write("\n}\n")

oFW.write("\ndictionaryReplacement\n{\n")
oFW.write("/* Sample entries\n")

oFW.write("\tboundary\n\t{\n")
oFW.write("\t\tFront\n\t\t{\n")
oFW.write("\t\t\ttype\twall;\n")
oFW.write("\t\t}\n")
oFW.write("\t}\n")

oFW.write("\tU\n\t{\n")
oFW.write("\t\tFront\n\t\t{\n")
oFW.write("\t\t\ttype\twall;\n")
oFW.write("\t\t}\n")
oFW.write("\t}\n")
oFW.write("}\n")
oFW.write("*/\n")

oFW.close()


# In[ ]:

## Create the files of the variables

variables = ["p","p_rgh","k","omega","epsilon","T","nut","U","alphat"]

for var in variables:
	oF = cwd + "/0/" + var
	oFW = open(oF,"w")

	for ii in xrange(0,len(header)):
		oFW.write(header[ii])

	if var == "U":
		oFW.write("\tclass\tvolVectorField;")
	else:
		oFW.write("\tclass\tvolScalarField;")
	oFW.write("\n\tobject\t"+var+";")
	oFW.write("\n}")
	
	if var == "p" or var == "p_rgh":
		oFW.write("\n\ndimensions\t[0 2 -2 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0;")
	elif var == "k":
		oFW.write("\n\ndimensions\t[0 2 -2 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "omega":
		oFW.write("\n\ndimensions\t[0 0 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "epsilon":
		oFW.write("\n\ndimensions\t[0 2 -3 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "T":
		oFW.write("\n\ndimensions\t[0 0 0 1 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t303.15;")
	elif var == "nut":	
		oFW.write("\n\ndimensions\t[0 2 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0;")
	elif var == "U":	
		oFW.write("\n\ndimensions\t[0 1 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t(1 0 0);")
	elif var == "alphat":	
		oFW.write("\n\ndimensions\t[0 2 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")

	oFW.write("\n\nboundaryField\n{\n")

	s1 = ["West","East","North"]
	for s in s1:
		oFW.write("\n\t"+s)
		oFW.write("\n\t{")
		if var == "p" or var == "p_rgh" or var == "T" or var == "alphat":
			oFW.write("\n\t\ttype\tzeroGradient;")
		elif var == "k" or var == "omega" or var == "epsilon":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform 0.1;")
		elif var == "nut":
			oFW.write("\n\t\ttype\tcalculated;")
			oFW.write("\n\t\tvalue\tuniform 0;")
		elif var == "U":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform (1 0 0);")	
		oFW.write("\n\t}")

	oFW.write("\n\tSouth")
	oFW.write("\n\t{")
	if var == "p" or var == "p_rgh":
		oFW.write("\n\t\ttype\tfixedValue;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var == "k" or var == "omega" or var == "epsilon":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "nut":
		oFW.write("\n\t\ttype\tcalculated;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var=="T":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "U":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "alphat":
		oFW.write("\n\t\ttype\tzeroGradient;")
	oFW.write("\n\t}")

	oFW.write("\n\tTop")
	oFW.write("\n\t{")
	if var == "p" or var == "T" or var == "alphat":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "p_rgh":
		oFW.write("\n\t\ttype\tfixedFluxPressure;")
		oFW.write("\n\t\tvalue\tinternalField;")
	elif var == "k" or var == "omega" or var == "epsilon":
		oFW.write("\n\t\ttype\tfixedValue;")
		oFW.write("\n\t\tvalue\tuniform 0.1;")
	elif var == "nut":
		oFW.write("\n\t\ttype\tcalculated;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var == "U":
		oFW.write("\n\t\ttype\tfixedValue;")
		oFW.write("\n\t\tvalue\tuniform (1 0 0);")	
	oFW.write("\n\t}")
	
	oFW.write("\n\tBottom")
	oFW.write("\n\t{")
	if var == "p":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "p_rgh":
		oFW.write("\n\t\ttype\tfixedFluxPressure;")
		oFW.write("\n\t\tvalue\tinternalField;")
	elif var == "k":
		oFW.write("\n\t\ttype\tkqRWallFunction;")
		oFW.write("\n\t\tvalue\tuniform 0.1;")
	elif var == "omega":
		oFW.write("\n\t\ttype\tomegaWallFunction;")
		oFW.write("\n\t\tvalue\tuniform 0.1;")
	elif var == "epsilon":
		oFW.write("\n\t\ttype\tepsilonWallFunction;")
		oFW.write("\n\t\tvalue\tuniform 0.1;")
	elif var == "nut":
		oFW.write("\n\t\ttype\tnutkWallFunction;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var=="T":
		#oFW.write("\n\t\ttype\tzeroGradient;")
		oFW.write("\n\t\ttype\tfixedValue;") #BB
		oFW.write("\n\t\tvalue\tuniform 0;") #BB
	elif var == "U":
		oFW.write("\n\t\ttype\tfixedValue;")
		oFW.write("\n\t\tvalue\tuniform (0 0 0);")
	elif var == "alphat":
		oFW.write("\n\t\ttype\talphatJayatillekeWallFunction;")
		oFW.write("\n\t\tPrt\t0.9;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	oFW.write("\n\t}")

	for jj in xrange(0,len(wSTL)):
		name = wSTL[jj]
		solid = wSolid[jj]
		oFW.write("\n\t"+name+"_"+solid)
		oFW.write("\n\t{")
		#if var == "p" or var == "T":
		if var == "p": #BB
			oFW.write("\n\t\ttype\tzeroGradient;")
		elif var == "p_rgh":
			oFW.write("\n\t\ttype\tfixedFluxPressure;")
			oFW.write("\n\t\tvalue\tinternalField;")
		elif var == "k":
			oFW.write("\n\t\ttype\tkqRWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;") 
		elif var == "omega":
			oFW.write("\n\t\ttype\tomegaWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;") 
		elif var == "epsilon":
			oFW.write("\n\t\ttype\tepsilonWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;")
		elif var == "nut":
			oFW.write("\n\t\ttype\tnutkWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0;")
		elif var == "U":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform (0 0 0);")	
		elif var == "alphat":
			oFW.write("\n\t\ttype\talphatJayatillekeWallFunction;")
			oFW.write("\n\t\tPrt\t0.9;")
			oFW.write("\n\t\tvalue\tuniform 0;")
 	        elif var == "T": #BB
			oFW.write("\n\t\ttype\ttimeVaryingHeatFlux;")
			oFW.write("\n\t\tgradient\tuniform 0;") 
			oFW.write("\n\t\theatSource\tflux;")
			oFW.write("\n\t\tq\tuniform 20;")
			oFW.write("\n\t\talphaEff\talphaEff;")
			oFW.write("\n\t\tvalue\t$internalField;")
		oFW.write("\n\t}")

	oFW.write("\n}")
	oFW.close()


# In[ ]:

## Create RASProperties and transportProperties files

# RASProperties
oF = cwd + "/constant/RASProperties"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tdictionary;")
oFW.write("\n\tlocation\tconstant;")
oFW.write("\n\tobject\tRASProperties;")
oFW.write("\n}\n")

oFW.write("\n\nRASModel\tkEpsilon;")
oFW.write("\n\nturbulence\ton;")
oFW.write("\n\nprintCoeffs\ton;")

oFW.close()

# transportProperties
oF = cwd + "/constant/transportProperties"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tdictionary;")
oFW.write("\n\tlocation\tconstant;")
oFW.write("\n\tobject\ttransportProperties;")
oFW.write("\n}\n")

oFW.write("\n\ntransportModel\tNewtonian;")
oFW.write("\n\nnu\tnu\t[0 2 -1 0 0 0 0]\t1.5E-05;")

oFW.write("\nTRef\tTRef\t[0 0 0 1 0 0 0]\t303.15;")
oFW.write("\nbeta\tbeta\t[0 0 0 -1 0 0 0]\t0.003298697;")
oFW.write("\nPr\tPr\t[0 0 0 0 0 0 0]\t0.9;")
oFW.write("\nPrt\tPrt\t[0 0 0 0 0 0 0]\t0.9;")
oFW.write("\n\nrhoCp0\t1005;")
for jj in xrange(0,len(wSTL)):
	oFW.write("\n"+wSTL[jj]+"_"+wSolid[jj]+".rhoCp0\t1005;")
oFW.close()

# Acceleration due to gravity g
oF = cwd + "/constant/g"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tuniformDimensionedVectorField;")
oFW.write("\n\tlocation\tconstant;")
oFW.write("\n\tobject\tg;")
oFW.write("\n}\n")

oFW.write("\n\ndimensions\t[0 1 -2 0 0 0 0];")
oFW.write("\n\nvalue\t(0 0 -9.81);")

oFW.close()


# In[ ]:

##Move STL files to constant/triSurface
oF = cwd + "/"
src = os.listdir(oF )
dst = cwd + "/constant/triSurface/"
for files in src:
	if files.endswith('.stl'):
		shutil.move(os.path.join(oF,files), os.path.join(dst,files))


# In[ ]:

os.mkdir("Misc")
os.rename("Input.VSGR","Misc/Input.VSGR")
os.rename("MasterSTLList","Misc/MasterSTLList")

try:
	os.mknod("case.foam")
except:
	oFW = open(cwd + "/case.foam","w")
	oFW.close()


# In[ ]:



