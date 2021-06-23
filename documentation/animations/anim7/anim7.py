# @brief Runs DFNMesh with rotating fracture. Trying to reproduce splinter bug in a small mesh
# @author github/PedroLima92

import subprocess
import math
import json
import numpy
import logging

#### Rotation matrices
# Qx = Rotation around x axis
def Qx(theta): 
    return numpy.matrix([[        1       ,        0       ,        0       ],
                         [        0       , math.cos(theta),-math.sin(theta)],
                         [        0       , math.sin(theta), math.cos(theta)]])

# Qy = Rotation around y axis
def Qy(theta): 
    return numpy.matrix([[ math.cos(theta),        0       , math.sin(theta)],
                         [        0       ,        1       ,        0       ],
                         [-math.sin(theta),        0       , math.cos(theta)]])

# Qz = Rotation around z axis
def Qz(theta): 
    return numpy.matrix([[ math.cos(theta),-math.sin(theta),        0       ],
                         [ math.sin(theta), math.cos(theta),        0       ],
                         [        0       ,        0       ,        1       ]])









# ../dfnMesh/<vtkfile.vtk>
vtkfile = "LOG/vtkmesh"
# where to save vtk files?
destination = "documentation/animations/anim7/"
# ../examples/<example>
inputpath = destination+"anim7.json"
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
data["comment"] = "An example to reproduce splinter formation on the surface of fractures. Likely caused by left over DFNFaces and DFNEdges that should've be removed once their polyhedron was refined"
data["Mesh"] = "examples/splinter.msh"
data["PreRefine"] = 1
data["TolDist"] = 10
data["TolAngle"] = 3.14159

nodesinitial = numpy.matrix.transpose(numpy.matrix([
    [-2.0, -2.0, -0.1],
    [ 2.0, -2.0, -0.1],
    [ 2.0,  2.0, -0.1],
    [-2.0,  2.0, -0.1]
]))
# nodes = numpy.matrix([
#     [-2.0, -2.0, -0.1],
#     [ 2.0, -2.0, -0.1],
#     [ 2.0,  2.0, -0.1],
#     [-2.0,  2.0, -0.1]
# ])


# Rotate fracture a few times
steps = 50
n_rotations_x = 1
n_rotations_y = 3*n_rotations_x
n_rotations_z = 3*n_rotations_y

step_angle_x = 2*math.pi*n_rotations_x/steps
step_angle_y = 2*math.pi*n_rotations_y/steps
step_angle_z = 2*math.pi*n_rotations_z/steps

# logging.debug("Step x = "+str(step_angle_x))
# logging.debug("Step y = "+str(step_angle_y))
# logging.debug("Step z = "+str(step_angle_z))

# print(nodes)

for i in range(42,43):
    theta_x = step_angle_x*i
    theta_y = step_angle_y*i
    theta_z = step_angle_z*i
    # logging.debug("theta x = "+str(theta_x))
    # logging.debug("theta y = "+str(theta_y))
    # logging.debug("theta z = "+str(theta_z))
    print(str(i)+" ==================")
    logging.debug("Rotation Matrix\n"+str(numpy.round(Qz(theta_z)*Qy(theta_y)*Qx(theta_x), decimals=7)))
    

    # Rotate fracture
    # nodes = Qx(theta_x)*Qy(theta_y)*Qz(theta_z)*nodesinitial
    nodes = (Qz(theta_z)*Qy(theta_y)*Qx(theta_x))*nodesinitial


    # Transpose node matrix to conform to DFNMesh syntax
    # nodelist = numpy.ndarray.tolist(numpy.round(numpy.matrix.transpose(nodes), decimals=4))
    nodelist = numpy.ndarray.tolist(numpy.matrix.transpose(nodes))
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
    dfnMesh = ["build/src/dfnTest"
              ,inputpath]
    subprocess.call(dfnMesh)
    # @todo exception handler to stop loop could go in here
    rename = ["cp"
              , vtkfile+".vtk"
              , destination+"vtkmesh"+"."+str(i)+".vtk"]
    subprocess.call(rename)

