# @brief Runs DFNMesh with variable tolerance and generate multiple vtk files to get a paraview animation
# @author github/PedroLima92
# @date june 2020

import subprocess

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/2020-06-22/"
# ../examples/<example>
example = destination+"variable-tol.txt"
# ../examples/<msh>
msh = "no-msh-file"

maxtoldist = 1.0
steps = 10
tol = maxtoldist/steps
steps += 1
# toldist = -tol

for i in range(steps):
    # run dfnMesh
    toldist = tol*i
    print(toldist)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,msh
              ,str(toldist)]
    subprocess.call(dfnMesh)
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)