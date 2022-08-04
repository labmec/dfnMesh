#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:03:31 2022

@author: Nathan
"""

# ----------- Imports -----------
import gmsh
import sys
import numpy as np

# ----------- Input -----------
# box domain
xmin = 0.
xmax = 1.
ymin = 0.
ymax = 2.25
zmin = 0.
zmax = 1.
elsize = 2.0 # not used if setting transfinite curve

# whre should it be cut
# ===> In x
xcuts = [0.05,0.167,0.334,0.5,0.667,0.833,0.95]

# ===> In y
ystart = 0.0; yend = 2.25
nycuts = 13
ycuts = np.array(range(nycuts-1)) + 1
ycuts = np.multiply(ycuts,(yend-ystart)/nycuts)
ycuts = np.insert(ycuts,0,0.05)
ycuts = np.insert(ycuts,ycuts.size,2.2)



# ===> In z
# NOTE: pay attention to boundary conditions when setting zcuts!!
zcuts = [0.05,0.167,0.334,0.5,0.667,0.833,0.95]

# Bounding box for inlet
xminIn = -0.01; xmaxIn = 1.01
yminIn = -0.01; ymaxIn = 0.01
zminIn = 1./3.-0.001; zmaxIn = 2./3.+0.001

# Bounding box for outlet
xminOut = -0.01; xmaxOut = 1.01
yminOut = 2.249; ymaxOut = 2.251
zminOut = -0.01; zmaxOut = 1./3.+0.001
zminOut2 = 2./3.-0.001; zmaxOut2 = 1.01

# ----------- Functions -----------
def create_plane(fixedPos: -1, fixedVal: float, elsize: float):
    assert(fixedPos == 0 or fixedPos == 1 or fixedPos == 2)

    vecpts = []    
    if fixedPos == 0:
        vecpts.append([fixedVal,ymin,zmin]); vecpts.append([fixedVal,ymax,zmin]); vecpts.append([fixedVal,ymax,zmax]); vecpts.append([fixedVal,ymin,zmax])
    elif fixedPos == 1:        
        vecpts.append([xmin,fixedVal,zmin]); vecpts.append([xmax,fixedVal,zmin]); vecpts.append([xmax,fixedVal,zmax]); vecpts.append([xmin,fixedVal,zmax])
    elif fixedPos == 2:       
        vecpts.append([xmin,ymin,fixedVal]); vecpts.append([xmax,ymin,fixedVal]); vecpts.append([xmax,ymax,fixedVal]); vecpts.append([xmin,ymax,fixedVal])
    
    pts = []

    for pt in vecpts:
        pts.append(gmsh.model.occ.add_point(pt[0], pt[1], pt[2], elsize))    

    l = []
    assert(len(pts) == 4)
    for i in range(len(pts)):
        p1 = i
        p2 = (i+1) % 4
        l.append(gmsh.model.occ.add_line(pts[p1], pts[p2]))
    lloop = gmsh.model.occ.add_curve_loop(l)
    surface = gmsh.model.occ.add_plane_surface([lloop])
    return surface

# --------------------------- Code -----------------------------
gmsh.initialize()

gmsh.model.add("fl_case3_coarse")

# Create initial cubic domain
gmsh.model.occ.addBox(xmin, ymin, zmin, xmax, ymax, zmax, 1)
gmsh.model.occ.synchronize()
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 2.0)

# Cut over predefined divisions x
fixedPosCut = 0
for c in xcuts:
    box = gmsh.model.getEntities(3)    
    surface = create_plane(fixedPosCut,c,elsize)
    gmsh.model.occ.fragment(box, [(2,surface)], removeObject=True, removeTool=True) # cuts volumes in box (Object) by surface (Tool)
    gmsh.model.occ.synchronize()

# Cut over predefined divisions y
fixedPosCut = 1
for c in ycuts:
    box = gmsh.model.getEntities(3)    
    surface = create_plane(fixedPosCut,c,elsize)
    gmsh.model.occ.fragment(box, [(2,surface)], removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

# Cut over predefined divisions z
fixedPosCut = 2
for c in zcuts:
    box = gmsh.model.getEntities(3)    
    surface = create_plane(fixedPosCut,c,elsize)
    gmsh.model.occ.fragment(box, [(2,surface)], removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

#  ========> Setting transfinite so it wont further divide each entity
for c in gmsh.model.getEntities(1):
    gmsh.model.mesh.setTransfiniteCurve(c[1], 2)
for s in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(s[1])
    gmsh.model.mesh.setRecombine(s[0], s[1])
    gmsh.model.mesh.setSmoothing(s[0], s[1], 2)
for v in gmsh.model.getEntities(3):
    gmsh.model.mesh.setTransfiniteVolume(v[1])    


# Creating Physical group for 3D volumes
dimToGet = 3
tagDom = 1
eps = 1e-4
vols = gmsh.model.getEntitiesInBoundingBox(xmin-eps,ymin-eps,zmin-eps,xmax+eps,ymax+eps,zmax+eps,dimToGet)
regions = [reg[1] for reg in vols]
nameDom = "k33"
gmsh.model.add_physical_group(dimToGet, regions, tagDom, nameDom)

# Creating Physical group for left inlet
dimToGet = 2
taginlet = 2
surfacesIn = gmsh.model.getEntitiesInBoundingBox(xminIn,yminIn,zminIn,xmaxIn,ymaxIn,zmaxIn,dimToGet)
regionsIn = [reg[1] for reg in surfacesIn]
nameIn = "inlet"
gmsh.model.add_physical_group(dimToGet, regionsIn, taginlet, nameIn)

# Creating Physical group for right outlets
dimToGet = 2
tagoutlet = 3
surfacesOut = gmsh.model.getEntitiesInBoundingBox(xminOut,yminOut,zminOut,xmaxOut,ymaxOut,zmaxOut,dimToGet)
surfaces2 = gmsh.model.getEntitiesInBoundingBox(xminOut,yminOut,zminOut2,xmaxOut,ymaxOut,zmaxOut2,dimToGet)
[surfacesOut.append(tup) for tup in surfaces2]
regionsOut = [reg[1] for reg in surfacesOut]
nameOut = "outlet"
gmsh.model.add_physical_group(dimToGet, regionsOut, tagoutlet, nameOut)

# No flux bcs on rest of surface
boundaries = gmsh.model.getBoundary(vols,True,False,False)
regionsNF = []
for reg in boundaries:
    thisTag = reg[1]
    isInletOrOutlet = thisTag in regionsIn or thisTag in regionsOut
    if not isInletOrOutlet:
        regionsNF.append(thisTag)
tagnoflux = 4
namenoflux = "noflux"
gmsh.model.add_physical_group(2, regionsNF, tagnoflux, namenoflux)

# Generate mesh:
gmsh.model.mesh.generate()
  
# Write mesh data:
gmsh.write("fl_case3_coarse.msh")

# Creates  graphical user interface
if 'close' not in sys.argv:
  gmsh.fltk.run()
  
# It finalize the Gmsh API
gmsh.finalize()
