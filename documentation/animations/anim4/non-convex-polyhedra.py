# @brief Runs DFNMesh with a mooving fracture in 3D which generates non-convex polyhedral regions
# @author github/PedroLima92
# @date july 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim4/"
# ../examples/<example>
example = destination+"anim4.txt"
# ../examples/<msh>
msh = "no-msh-file"

preamble = (
"""Domain
4.0 6.0 4.0

Mesh
EHexahedral
1 1 1
"""
)
x = "\n-1.00  5.00  5.00 -1.00"
z = "\n-1.00 -1.00  5.00  5.00"

toldist = 0.6
# toldist = 0.00006
preamble += "\nTolDist  " + str(toldist) + "\n"
step = 0.05
# steps = 0
yinitial = -1.0
yfinal = 6.5
steps = int((yfinal-yinitial)/step)
steps += 1

# for i in range(steps):
for i in range(38):
# for i in range(32):
# for i in range(14,15):
    print("\nstep ",i,"\n")
    y0 = yinitial+step*i
    f = open(example,"w+")
    f.write(preamble)
    f.write("\nFracture 0")
    f.write("\nLimit Etruncated")
    f.write(x)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0-0.7,y0,y0+0.7,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0-0.7,y0,y0+0.7,y0))
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