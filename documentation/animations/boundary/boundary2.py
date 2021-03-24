# @brief Example to test the fracture boundary recovery
# @author github/PedroLima92
# @date november 2020

import subprocess
# import math
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import animtools









# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/boundary/ex2/"
# ../examples/<example>
example = destination+"boundary2.txt"
# ../examples/<msh>
msh = "no-msh-file"






# grid origin
origin = [ 0.0,  0.0,  0.0]
# domain dimensions
domain = [ 4.0,  4.0,  4.0]
# mesh partitions
eltype = "EHexahedral"
partitions = [2, 2, 2]

# preamble
preamble = ""
if(msh == "no-msh-file"):
    preamble += "Origin\n"+str(origin[0])+' '+str(origin[1])+' '+str(origin[2])+' '+"\n\n"
    preamble += "Domain\n"+str(domain[0])+' '+str(domain[1])+' '+str(domain[2])+' '+"\n\n"
    preamble += "Mesh\n"+eltype+"\n"+str(partitions[0])+' '+str(partitions[1])+' '+str(partitions[2])+' '+"\n\n"
else:
    preamble += "Mesh\n\""+msh+"\""

preamble += "\n\n"

# Fractures
fracture0 = "Fracture 0\n"
fracture0 += "Limit Erecovered\n"
# p0 = [[  0.1,  0.8,  0.8,  0.1],
#       [  7.0, -0.1, -0.1,  7.0],
#       [ -1.0, -1.0,  1.0,  1.0]]
# for i in range(3):
#     for j in range(len(p0[0])):
#         if(p0[i][j] >= 0.0): fracture0 += ' '
#         fracture0 += str(p0[i][j])+' '
#     fracture0 += "\n"

fracture1 = "Fracture 1\n"
fracture1 += "Limit Erecovered\n"
# p0 = [[  7.0, -0.1, -0.1,  7.0],
#       [  0.1,  0.8,  0.8,  0.1],
#       [ -1.0, -1.0,  1.0,  1.0]]
# for i in range(3):
#     for j in range(len(p0[0])):
#         if(p0[i][j] >= 0.0): fracture1 += ' '
#         fracture1 += str(p0[i][j])+' '
#     fracture1 += "\n"


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
    # f.write("\n\nFracture 0 4")
    # f.write(x0)
    if y0 < 0 :
        y = ("\n%.2f %.2f %.2f %.2f" % (y0,y0,y0,y0))
    else:
        y = ("\n %.2f  %.2f  %.2f  %.2f" % (y0,y0,y0,y0))
    # f.write(y)
    # f.write(z0)
    f.write(fracture0)
    f.write(y)
    f.write(x0)
    f.write(z0)
    f.write("\n\n")
    f.write(fracture1)
    f.write(x1)
    f.write(y)
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