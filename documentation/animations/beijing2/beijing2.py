# @brief Example number 2 for LSEC Seminar 2020 (Beijing)
# @author github/PedroLima92
# @date october 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/beijing2/"
# ../examples/<example>
example = destination+"beijing2.txt"
# ../examples/<msh>
msh = "examples/tetrahedra2.msh"




firstFracture = (
"""\n\n
Fracture 1 4
 0.50  0.50  0.50  0.50
-1.50  1.50  1.50 -1.50
-1.50 -1.50  1.50  1.50
""")





preamble = ("NumberOfFractures 1")
# """Domain
# 4.0 4.0 4.0

# Mesh
# EQuadrilateral
# 2 2 2

# NumberOfFractures 2

# """
# )

x = "\n-2.00  2.00  2.00 -2.00"
z = "\n-2.00 -2.00  2.00  2.00"

toldist = 0.2
tolangle = 0.5
step = 0.07
# steps = 0
starty0 = -1.0
yfinal = 0.9999
steps = int((yfinal-starty0)/step)
steps += 1

for i in range(steps):
# for i in range(3,4):
    print("step ",i)
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
    # f.write("\n\nFracture 1 4")
    # f.write(y)
    # f.write(x)
    # f.write(z)
    f.close()
    # run dfnMesh
    print(y0)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,msh
              ,str(toldist)
              ,str(tolangle)]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)


# for i in range(steps):