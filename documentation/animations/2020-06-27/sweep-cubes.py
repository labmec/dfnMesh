# @brief Runs DFNMesh with a mooving fracture in 3D which causes some 2D elements to be incorporated into the fracture
# @author github/PedroLima92
# @date june 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/2020-06-27/"
# ../examples/<example>
example = destination+"sweep-cubes.txt"
# ../examples/<msh>
msh = "no-msh-file"

preamble = (
"""Domain
4.0 4.0 4.0

Mesh
EQuadrilateral
2 2 2

NumberOfFractures 1

Fracture 0 4"""
)
x = "\n-1.00  5.00  5.00 -1.00"
z = "\n-1.00 -1.00  5.00  5.00"

toldist = 0.6
step = 0.05
# steps = 0
y0 = -1.0
yfinal = 5.0
steps = int((yfinal-y0)/step)
steps += 1

for i in range(steps):
    f = open(example,"w+")
    f.write(preamble)
    f.write(x)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0,y0,y0,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0,y0,y0,y0))
    f.write(y)
    f.write(z)
    f.close()
    # run dfnMesh
    print(y0)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,msh
              ,str(toldist)]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)
    y0 += step


# for i in range(steps):