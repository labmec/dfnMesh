

<p align="center"> 
 <h1> DFNMesh </h1>
</p>

<p align="justify">
The main goal of this project is to provide an intersection searching-and-processing tool that is able to insert polygons into a tridimensional coarse mesh and refine it to conform to these polygons. Ultimately describing a fractured medium, as a two-level discretisation mesh, for Discrete Fracture Network (DFN) flow simulation methods available in <a href="https://github.com/labmec/neopz">NeoPZ</a>.
 </p>

<p align="center">
<img src="./documentation/img/octagon-macros.png" title="Fracture planes in coarse mesh" width="30%" height="30%"/>
<img src="./documentation/img/octagon-frac.png" title="Fracture surfaces" width="25%" height="25%" /> 
<img src="./documentation/img/octagon-vol.png" title="Submesh" width="30%" height="30%" />
</p>

## Dependencies
This project is highly dependent in 2 open source libs:
### [NeoPZ](https://github.com/labmec/neopz "NeoPZ repository")
### [GMsh](https://gitlab.onelab.info/gmsh/gmsh "GMsh repository")

(see <code>main()</code> in each branch for its assumed version of each library)

## Cool graphics
<p>
<img src="./documentation/img/colored-refpatterns.png" title="Refinement patterns conforming to fracture surface" width="40%" height="40%"/>
</p>

