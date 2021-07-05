# @brief Runs DFNMesh with translating fracture over a semi-spherical domain
# @author github/PedroLima92

import subprocess
import json
import numpy
import logging









# ../dfnMesh/<vtkfile.vtk>
vtkfile = "LOG/vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim8/"
# ../examples/<example>
inputpath = destination+"anim8.json"
# ../examples/<msh>
# msh = "examples/splinter.msh"
msh = destination+"/testmesh.msh"

# Config log file
logging.basicConfig(filename=destination+"logfile.log",
                    level=logging.DEBUG,
                    filemode='w',
                    format="%(levelname)s - %(message)s"
)



# Setup initial json data for input file
data = {}
data["$schema"] = "../../../examples/dfn_schema.json"
data["comment"] = "Animation for UFPE seminar"
data["Mesh"] = "examples/sphere.msh"
data["PreRefine"] = 0
data["TolDist"] = 10
data["TolAngle"] = 3.14159

nodesinitial = numpy.matrix([
    [-1.05,-1.05,-0.05],
    [-1.05, 1.05,-0.05],
    [-1.05, 1.05, 1.05],
    [-1.05,-1.05, 1.05]
])
nodesfinal = numpy.matrix([
    [ 1.05,-1.05,-0.05],
    [ 1.05, 1.05,-0.05],
    [ 1.05, 1.05, 1.05],
    [ 1.05,-1.05, 1.05]
])


# Translate fracture along x
steps = 50
dx = (nodesfinal[0,0] - nodesinitial[0,0])/steps
translation = numpy.full(nodesinitial.shape, 0.0)
for i in range(nodesinitial.shape[0]):
    translation[i,0] = dx

logging.info("Step translation vector:\n" + str(translation))


for i in range(steps):
    print(str(i)+" ==================")
    

    # Translate fracture
    nodes = nodesinitial + translation*i
    logging.debug("Nodes, step "+str(i)+'\n'+str(numpy.round(nodes, decimals=7)))


    # Add fracture to dictionary
    nodelist = numpy.ndarray.tolist(nodes)
    data["Fractures"] = []
    data["Fractures"].append({
        "Index": 0,
        "Limit": "Eextended",
        "MaterialID": 10,
        "Nodes": nodelist
    })

    # Write json to input file
    with open(inputpath, 'w') as inputfile:
        json.dump(data, inputfile, indent=4)

    # run dfnMesh
    dfnMesh = ["build/src/animations"
              ,inputpath]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+"vtkmesh"+"."+str(i)+".vtk"]
    subprocess.call(rename)

