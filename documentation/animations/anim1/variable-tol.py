# @brief Runs DFNMesh with variable tolerance and generate multiple vtk files to get a paraview animation
# @author github/PedroLima92
# @date june 2020

import subprocess

def MatrixToString(fmatrix):
    rows = len(fmatrix)
    cols = len(fmatrix[0])
    fstring = ""
    for i in range(rows):
        for j in range(cols):
            if(p0[i][j] >= 0.0): fstring += ' '
            fstring += str(p0[i][j])+' '
        fstring += "\n"
    print(fstring)
    return fstring


# ../dfnMesh/<vtkfile.vtk>
vtkfile = "vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim1/"
# ../examples/<example>
example = destination+"anim1.txt"
# ../examples/<msh>
msh = "no-msh-file"






# grid origin
origin = [-3.0,  0.0,  0.0]
# domain dimensions
domain = [ 6.0,  8.0,  0.0]
# mesh partitions
eltype = "EQuadrilateral"
partitions = [3, 4, 0]

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
fracture0 += "Limit Etruncated\n"
p0 = [[  0.1,  0.8,  0.8,  0.1],
      [  7.0, -0.1, -0.1,  7.0],
      [ -1.0, -1.0,  1.0,  1.0]]
fracture0 += MatrixToString(p0)

maxtoldist = 1.0
steps = 10
tol = maxtoldist/steps
steps += 1

for i in range(steps):
    toldist = tol*i
    myfile = open(example,"w+")
    myfile.write(preamble)
    myfile.write("TolDist "+str(toldist)+"\n\n")
    myfile.write(fracture0+"\n")
    myfile.close()
    # run dfnMesh
    dfnMesh = ["build/src/dfnTest"
              ,"-f",example
              ,"-m",msh
              ,"-td",str(toldist)]
    subprocess.call(dfnMesh)
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+vtkfile+"."+str(i)+".vtk"]
    subprocess.call(rename)