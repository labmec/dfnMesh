<h1 align="center">DFNMesh</h1>

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
<img src="https://s3.us-west-2.amazonaws.com/secure.notion-static.com/e99a63e7-1b61-4334-b19a-52d0d5120413/2020-06-26_18-11-02.gif?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAT73L2G45O3KS52Y5%2F20200627%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20200627T194433Z&X-Amz-Expires=86400&X-Amz-Signature=f1dc4ec034c575a9a5f26a9774d01f4894b34f48618e0bdef218a8b65a20658c&X-Amz-SignedHeaders=host" title="Coalescing fracture that is close to nodes/edges" width="40%" height="40%"/>
</p>
