# @brief Runs DFNMesh with rotating fracture (around z-axis) in 3D over a mesh of tetrahedra
# @author github/PedroLima92
# @date july 2020

import subprocess
import math

# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/2020-07-22/"
# ../examples/<example>
example = destination+"tetrahedra1.txt"
# ../examples/<msh>
msh = "examples/tetrahedra1.msh"

preamble = (
"""NumberOfFractures 1

Fracture 0 4"""
)
z = "\n-1.50 -1.50  1.50  1.50"

toldist = [0.01, 0.333]
# toldist = [0.01, 0.05]
tolangle = [0.01, (19)*math.pi/180]
step = 0.5
# steps = 0
steps = int(180/step)
step = step*math.pi/180
steps += 1
x0 = -0.3
y0 = -0.3
r = 2.0
theta = 0.0

for i in range(steps):
# for i in range(67,68):
    # theta = (i-i%2)*step
    theta = i*step
    f = open(example,"w+")
    f.write(preamble)
    x1 = x0-r*math.cos(theta)
    y1 = y0-r*math.sin(theta)
    x2 = x0+r*math.cos(theta)
    y2 = y0+r*math.sin(theta)
    x = ("\n %.2f %.2f %.2f %.2f" % (x1,x2,x2,x1))
    y = ("\n %.2f %.2f %.2f %.2f" % (y1,y2,y2,y1))
    f.write(x)
    f.write(y)
    f.write(z)
    f.close()
    # run dfnMesh
    print(toldist)
    dfnMesh = ["build/src/dfnTest"
              ,example
              ,msh
              ,str(toldist[1])
              ,str(tolangle[1])]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)


# for i in range(steps):