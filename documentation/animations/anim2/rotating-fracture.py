# @brief Runs DFNMesh with rotating fracture in 2D and generate multiple vtk files to get a paraview animation
# @author github/PedroLima92
# @date june 2020

import subprocess
import math




def MatrixToString(fmatrix):
    rows = len(fmatrix)
    cols = len(fmatrix[0])
    fstring = ""
    for i in range(rows):
        for j in range(cols):
            if(fmatrix[i][j] >= 0.0): fstring += ' '
            fstring += str(fmatrix[i][j])+' '
        fstring += "\n"
    return fstring






# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim2/"
# ../examples/<example>
example = destination+"rotating-fracture.txt"
# ../examples/<msh>
msh = "no-msh-file"





# grid origin
origin = [0.0,  0.0,  0.0]
# domain dimensions
domain = [ 8.0,  8.0,  0.0]
# mesh partitions
eltype = "EQuadrilateral"
partitions = [4, 4, 0]

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

x0 = [0.,0.,0.]
x0[0] = 0.5*(domain[0]-origin[0])
x0[1] = 0.5*(domain[1]-origin[1])
x0[2] = 0.5*(domain[2]-origin[2])

fracture0 = "Fracture 0\n"
fracture0 += "Limit Etruncated\n"
preamble += fracture0
# p0 = [[  0.1,  0.8,  0.8,  0.1],
#       [  7.0, -0.1, -0.1,  7.0],
#       [ -1.0,  1.0,  1.0, -1.0]]
# for i in range(3):
#     for j in range(len(p0[0])):
#         if(p0[i][j] >= 0.0): fracture0 += ' '
#         fracture0 += str(p0[i][j])+' '
#     fracture0 += "\n"





z = "\n-1.00  1.00  1.00 -1.00"

toldist = 0.5
step = 0.5
steps = int(180/step)
step = step*math.pi/180
steps += 1
x0 = 4.0
y0 = 4.0
r = 5.66
theta0 = 0

for i in range(steps):
    theta = theta0 + i*step
    f = open(example,"w")
    f.write(preamble)
    
    x1 = x0+r*math.cos(theta)
    y1 = y0+r*math.sin(theta)
    x2 = x0-r*math.cos(theta)
    y2 = y0-r*math.sin(theta)
    x = ("\n% .2f % .2f % .2f % .2f" % (x1,x1,x2,x2))
    y = ("\n% .2f % .2f % .2f % .2f" % (y1,y1,y2,y2))
    f.write(x)
    f.write(y)
    f.write(z)
    f.close()
    # run dfnMesh
    print(toldist)
    dfnMesh = ["build/src/dfnTest"
              ,"-f",example
              ,"-m",msh
              ,"-td",str(toldist)]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)


# for i in range(steps):