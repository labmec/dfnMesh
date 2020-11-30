# @brief Runs DFNMesh with a mooving fracture in 3D which causes some 2D elements to be incorporated into the fracture
# @author github/PedroLima92
# @date june 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim3/"
# ../examples/<example>
example = destination+"anim3.txt"
# ../examples/<msh>
msh = "no-msh-file"

preamble = (
"""Domain
4.0 4.0 4.0

Mesh
EHexahedral
2 2 2

""")

toldist = 0.8
preamble += "TolDist  "+str(toldist)+"\n"



x = "\n-1.00  5.00  5.00 -1.00"
z = "\n-1.00 -1.00  5.00  5.00"
step = 0.07
# steps = 0
y0 = 1.85
yfinal = 5.0
steps = int((yfinal-y0)/step)
steps += 1

# for i in range(17,18):
for i in range(steps):
    y0 += step
    f = open(example,"w+")
    f.write(preamble)
    f.write("\nFracture 0")
    f.write("\nLimit Etruncated")
    f.write(x)
    y = ("\n% .2f % .2f % .2f % .2f" % (y0,y0,y0,y0))
    f.write(y)
    f.write(z)
    f.close()
    # run dfnMesh
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