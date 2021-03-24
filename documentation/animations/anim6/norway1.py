# @brief Multiple fractures in a single cell to demonstrate conservation of original faces and nodes
# @author github/PedroLima92
# @date november 2020

import subprocess
# import math
# import os,sys,inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0,parentdir) 

# import animtools


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
destination = "documentation/animations/anim6/out1/"
# ../examples/<example>
example = destination+"norway1.txt"
# ../examples/<msh>
msh = "no-msh-file"






# grid origin
origin = [-1.0, -1.0, -1.0]
# domain dimensions
domain = [ 2.0,  2.0,  2.0]
# mesh partitions
eltype = "EHexahedral"
partitions = [1, 1, 1]

# preamble
preamble = ""
if(msh == "no-msh-file"):
    preamble += "Origin\n"+str(origin[0])+' '+str(origin[1])+' '+str(origin[2])+' '+"\n\n"
    preamble += "Domain\n"+str(domain[0])+' '+str(domain[1])+' '+str(domain[2])+' '+"\n\n"
    preamble += "Mesh\n"+eltype+"\n"+str(partitions[0])+' '+str(partitions[1])+' '+str(partitions[2])+' '+"\n\n"
else:
    preamble += "Mesh\n\""+msh+"\""

# preamble += "\n"

# Fractures
fracture0 = "Fracture 0\n"
fracture0 += "Limit Etruncated\n"
p0 = [[ -1.5,  1.5,  1.5, -1.5],
      [ -1.5, -1.5,  1.5,  1.5],
      [  0.0,  0.0,  0.0,  0.0]]
fracture0 += MatrixToString(p0)

fracture1 = "Fracture 1\n"
fracture1 += "Limit Etruncated\n"
p1 = [[  0.0,  0.0,  0.0,  0.0],
      [ -1.5,  1.5,  1.5, -1.5],
      [ -1.5, -1.5,  1.5,  1.5]]
fracture1 += MatrixToString(p1)

fracture2 = "Fracture 2\n"
fracture2 += "Limit Etruncated\n"
p2 = [[ -1.5,  1.5,  1.5, -1.5],
      [  0.0,  0.0,  0.0,  0.0],
      [ -1.5, -1.5,  1.5,  1.5]]
fracture2 += MatrixToString(p2)

fracture3 = "Fracture 3\n"
fracture3 += "Limit Etruncated\n"
p3 = [[ -1.5,  1.5,  1.5, -1.5],
      [ -1.5, -1.5,  1.5,  1.5],
      [  1.5,  1.5, -1.5, -1.5]]
fracture3 += MatrixToString(p3)

fracture4 = "Fracture 4\n"
fracture4 += "Limit Erecovered\n"
p4 = [[ 0.5,  1.6,  1.6,  0.5],
      [ 0.5,  1.6,  1.6,  0.5],
      [ 0.5,  0.5,  1.6,  1.6]]
fracture4 += MatrixToString(p4)


fractures = [fracture0,fracture1,fracture2,fracture3,fracture4]
nfrac = len(fractures)

toldist = 0.2

file = open(example,"w")
file.write(preamble)
file.write("TolDist  " + str(toldist) + "\n")
# file.write("TolAngle " + str(tolangle) + "\n")
file.close()

for i in range(-1, nfrac + 1):
    file = open(example,"a")
    if(i < nfrac and i > -1):
        file.write("\n" + fractures[i])
    file.close()


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
              , destination+vtkfile+"."+str(i+1)+".vtk"]
    subprocess.call(rename)

