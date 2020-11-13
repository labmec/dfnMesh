# @brief Example number 2 for cilamce2020 paper
# @author github/PedroLima92
# @date august 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/cilamce/ex2/"
# ../examples/<example>
example = destination+"cilamce2.txt"
# ../examples/<msh>
msh = "no-msh-file"




firstFracture = (
"""\n\n
Fracture 1 4
 3.00  3.00  3.00  3.00
-1.00  5.00  5.00 -1.00
-1.50 -1.50  5.50  5.50
""")





preamble = (
"""Domain
4.0 4.0 4.0

Mesh
EHexahedral
2 2 2

NumberOfFractures 2

"""
)
x = "\n-1.00  5.00  5.00 -1.00"
z = "\n-1.00 -1.00  5.00  5.00"

toldist = 0.4
step = 0.07
# steps = 0
starty0 = -0.1
yfinal = 5.0
steps = int((yfinal-starty0)/step)
steps += 1

for i in range(steps):
# for i in range(17,18):
    y0 = starty0 + i*step
    f = open(example,"w+")
    f.write(preamble)
    f.write("\n\nFracture 0 4")
    f.write(x)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0,y0,y0,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0,y0,y0,y0))
    f.write(y)
    f.write(z)
    f.write("\n\nFracture 1 4")
    f.write(y)
    f.write(x)
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