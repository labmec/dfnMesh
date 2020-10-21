# @brief Runs DFNMesh with a mooving fracture in 3D which generates non-convex polyhedral regions
# @author github/PedroLima92
# @date july 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/2020-07-20/"
# ../examples/<example>
example = destination+"non-convex-polyhedra.txt"
# ../examples/<msh>
msh = "no-msh-file"

preamble = (
"""Domain
4.0 6.0 4.0

Mesh
EQuadrilateral
2 2 2

NumberOfFractures 1

Fracture 0 4"""
)
x = "\n-1.00  5.00  5.00 -1.00"
z = "\n-1.00 -1.00  5.00  5.00"

toldist = 0.6
# toldist = 0.00006
step = 0.05
# steps = 0
yinitial = -1.0
yfinal = 6.5
steps = int((yfinal-yinitial)/step)
steps += 1

for i in range(steps):
# for i in range(101,102):
    y0 = yinitial+step*i
    f = open(example,"w+")
    f.write(preamble)
    f.write(x)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0-0.7,y0,y0+0.7,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0-0.7,y0,y0+0.7,y0))
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


# for i in range(steps):