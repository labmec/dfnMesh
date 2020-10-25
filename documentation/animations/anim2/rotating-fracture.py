# @brief Runs DFNMesh with rotating fracture in 2D and generate multiple vtk files to get a paraview animation
# @author github/PedroLima92
# @date june 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/2020-06-23/"
# ../examples/<example>
example = destination+"rotating-fracture.txt"
# ../examples/<msh>
msh = "no-msh-file"

preamble = (
"""Domain
8.0 8.0 0.0

Mesh
EQuadrilateral
4 4 0

NumberOfFractures 1

Fracture 0 4"""
)
z = "\n-1.00  1.00  1.00 -1.00"

toldist = 0.5
step = 0.5
# steps = 0
steps = int(180/step)
step = step*math.pi/180
steps += 1
x0 = 4.0
y0 = 4.0
r = 5.66
theta0 = 0

for i in range(steps):
    theta = theta0 + i*step
    f = open(example,"w+")
    f.write(preamble)
    x1 = x0+r*math.cos(theta)
    y1 = y0+r*math.sin(theta)
    x2 = x0-r*math.cos(theta)
    y2 = y0-r*math.sin(theta)
    x = ("\n %.2f %.2f %.2f %.2f" % (x1,x1,x2,x2))
    y = ("\n %.2f %.2f %.2f %.2f" % (y1,y1,y2,y2))
    f.write(x)
    f.write(y)
    f.write(z)
    f.close()
    # run dfnMesh
    print(toldist)
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