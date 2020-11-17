# @brief Example to test the fracture boundary recovery
# @author github/PedroLima92
# @date november 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/boundary/ex1/"
# ../examples/<example>
example = destination+"boundary1.txt"
# ../examples/<msh>
msh = "no-msh-file"




preamble = (
"""Domain
4.0 4.0 4.0

Mesh
EHexahedral
2 2 2

"""
)
x0 = "\n-1.50  5.00  5.00 -1.50"
z0 = "\n-0.50 -0.50  5.00  5.00"
x1 = "\n 1.50  5.00  5.00  1.50"
z1 = "\n 1.50  1.50  5.00  5.00"

toldist = 0.4
step = 0.07
# steps = 0
starty0 = -0.1
yfinal = 5.0
steps = int((yfinal-starty0)/step)
steps += 1

for i in range(steps):
# for i in range(12,13):
    y0 = starty0 + i*step
    f = open(example,"w+")
    f.write(preamble)
    f.write("\n\nFracture 0 4")
    f.write(x0)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0,y0,y0,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0,y0,y0,y0))
    f.write(y)
    f.write(z0)
    f.write("\n\nFracture 1 4")
    f.write(y)
    f.write(x1)
    f.write(z1)
    f.close()
    # run dfnMesh
    print(y0)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,"-m"
              ,msh
              ,"-td"
              ,str(toldist)]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)


# for i in range(steps):