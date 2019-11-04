# Ongoing

- [ ] Geometrical center -> Centroid
- [ ] Insert second plane
  - [x] Implement triangle splitting cases
  - [x] Set a second plane in main
  - [x] Call splitribs and splitfaces
  - [x] Set RefPatterns
    - [x] Refpattern for ribs shouldn't be uniform
  - [x] Use PZ's data structure to specify level of element splitting
- [ ] Tetrahedralization on volumetric mesh (gmesh should do it) 
  - [x] Search which volume encloses each face
  - [x] Search through 2nd (& 3rd, 4th...) level neighbours until a child of an MHM element is found
  - [x] Inform to volumes what faces lie inside it 
  - [ ] **Write geometry to GMSH** 
    - [x] Write points 
    - [x] Write Lines 
    - [x] Define physical groups of lines 
    - [x] Write faces 
    - [x] Define physical groups of faces 
    - [x] Test rules for proper definition of surface loops 
    - [x] Write surface loops 
    - [x] Write volumes 
    - [x] Define physical groups of volumes 
    - [x] Rewrite using gmsh API
    - [x] One volume at a time (to allow for refinement patterns)
    - [ ] Check if child is of lower dimension than father
    - [ ] When sending stuff to GMsh, there might not be any intact volume, in which case, they must be ignored.
# ToDo
- [X] Swap pzgeoel::NNodes() for pzgeoel::NCornerNodes where necessary
- [x] Mesh fracplane 
  - [X] ~~*Ear-clipping-like algorithm*~~ [2019-09-30] 
  - [x] Legalize triangles 
  - [x] GMsh is doing it perfectly
- [ ] Face splitting 
  - [x] Store status of each rib and node
  - [x] Identify case 
  - [x] Store all intersections for each face 
  - [x] Father elements 
  - [ ] Quadrilateral splitting 
    - [x] Case 1 
    - [x] Case 2 
    - [ ] Case 3 
    - [ ] Case 4
    - [x] Case 5 
    - [ ] Case 6 
    - [ ] Case 7 
    - [ ] Case 8 (No ribs 1) 
    - [ ] Case 9 (No ribs 2) 
  - [ ] Triangle splitting 
    - [x] Case 10 
    - [ ] Case 11 
    - [x] Case 12 
    - [ ] Case 13 
    - [ ] Case 14 
    - [ ] Case 15 (No ribs 1) 
    - [ ] Case 16 (No ribs 2) 
- [ ] Tetrahedralization on volumetric mesh (gmesh should do it) 
  - [x] Search which volume encloses each face
  - [x] Inform to volumes what faces lie inside it 
  - [x] Write geometry to GMSH 
    - [x] Write points 
    - [x] Write Lines 
    - [x] Define physical groups of lines 
    - [x] Write faces 
    - [x] Define physical groups of faces 
    - [x] Test rules for proper definition of surface loops 
    - [x] Write surface loops 
    - [x] Write volumes 
    - [x] Define physical groups of volumes 
    - [x] Rewrite using gmsh API
- [ ] Check if rib cut happens too close to vertex 
- [ ] Check if face intersection happens too close to rib 
- [ ] Search for critical cases 
  - [ ] Fractures that are almost parallel to mesh plane 
  - [ ] Fracture that overlaps mesh plane should not divide anything 
- [ ] Id volumes that contain fracture corners 
- [ ] Check if point is too close to vertex -> rib -> face 
- [ ] Change comparisons from a == b to $|a-b|<\varepsilon$


# Done and Documented
- [x] Rib class
- [x] Face class
- [x] Fracture Plane class
  - [x] Check data consistency
  - [x] Compute axis
  - [x] Compute area
- [x] Fracture Mesh class
  - [x] Map of ribs
  - [x] Map of Mid-Faces
  - [x] Map of End-Faces
- [x] Create skeleton mesh
  - [x] Check if element already has skeleton 
- [x] Check if point is above or bellow plane
- [x] Check if rib is cut by infinite plane
- [x] Calculate intersection with inifinite plane
- [x] Check if intersection is within plane boundaries
  - [x] Compute area of triangles from every corner of plane to point
  - [x] Compare area of triangles to area of plane (to a tolerance)
- [x] Find faces that are cut
  - [x] Iterate over faces
  - [x] Iterate over its 1D neighbours through 1D sides
  - [x] Store in mid-faces if has 2 ribs cut
  - [x] Store in end-faces if has 1 rib cut
- [x] Find intersection points in end-faces
  - [x] figure out 'prettiest' conversion from TPZGeoEl to DFNFracPlane
  - [x] create EPoint for intersection points 
- [x] Connect intersection points to define fracture surface outline
  - [x] Connect them at the time an intersected face is found and pair of points are created
- [x] Split fracture front (fracture surface edges)
  - [x] iterate over end-faces
  - [x] find which edge cuts them
  - [x] sort them
  - [x] connect points accordingly






# Archived
Phil's suggestions
- [ ] Search for volumes that contain corner nodes
  - [ ] Ribs that have that node, are contained by the same volume
  - [ ] Faces that have that node, are contained by the same volume
- [ ] From the fracture outline, contours may be generated and these have an easily found containing volume. All elements in that contour will inevitably belong to the same volume.
- [ ] If eldest ancestor has father ID of an element of higher dimension than itself, then that father ID belongs to a MHM element and that MHM element is the enclosing volume (or enclosing area).
- [ ] If we use fatherID to define enclosing volume, would we hurt the robustness of PZ's MHM methods?

Fran's suggestions
- [x] Search through 2nd (& 3rd, 4th...) level neighbours until a child of an MHM element is found
- [ ] Maybe 2D elements could carry the index of their enclosing volumes, this could be a useful optimization
    if neighbour has enclosing volume
        test for that enclosing volume
        if not found
            test for all 3D neighbours