# @brief Runs DFNMesh with variable tolerance and generate multiple vtk files to get a paraview animation
# @author Pedro Lima
# @date june 2020

import sh

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# ../examples/<example>
example = "variable-tol.txt"
# ../examples/<msh>
msh = "skip"

maxtoldist = 1.0
steps = 5
tol = maxtoldist/steps


for i in range(steps):
    # run dfnMesh
    rename = ["cp"
              , vtkfile+".vtk"
              , vtkfile+"."+str(i)+".vtk"]
    tol += tol