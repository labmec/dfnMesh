# @brief Runs DFNMesh with variable tolerance and generate multiple vtk files to get a paraview animation
# @author github/PedroLima92
# @date june 2020

import subprocess

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# ../examples/<example>
example = "examples/variable-tol.txt"
# ../examples/<msh>
msh = "no-msh-file"

maxtoldist = 1.5
steps = 10
tol = maxtoldist/steps
toldist = -tol

for i in range(steps):
    # run dfnMesh
    toldist += tol
    print(toldist)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,msh
              ,str(toldist)]
    subprocess.call(dfnMesh)
    rename = ["cp"
              , vtkfile+".vtk"
              , vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)