# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
vtkmesh = FindSource('vtkmesh.*')

# ELEMENT RANGE
elementRange = Threshold(Input=vtkmesh)
elementRange.Scalars = ['CELLS', 'elIndex']
elementRange.ThresholdRange = [0.0, 10000.0]
RenameSource('ElementRange', elementRange)



# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# get layout
layout1 = GetLayout()


# 3D FULL
dim3Dfull = Threshold(Input=elementRange)
dim3Dfull.Scalars = ['CELLS', 'Dimension']
dim3Dfull.ThresholdRange = [2.0, 3.0]
RenameSource('3D', dim3Dfull)

# 1D FULL
dim1Dfull = Threshold(Input=elementRange)
dim1Dfull.Scalars = ['CELLS', 'Dimension']
dim1Dfull.ThresholdRange = [1.0, 1.0]
RenameSource('1D', dim1Dfull)

# COARSE MESH
coarseMesh = Threshold(Input=dim3Dfull)
coarseMesh.Scalars = ['CELLS', 'elIndex']
coarseMesh.ThresholdRange = [0.0, 16.0]
RenameSource('CoarseMesh', coarseMesh)
# show data in view
coarseMeshDisplay = Show(coarseMesh, renderView1, 'UnstructuredGridRepresentation')
coarseMeshDisplay.SetRepresentationType('Wireframe')
coarseMeshDisplay.Opacity = 0.3

# iPOINTS
iPoints = Threshold(Input=dim1Dfull)
iPoints.Scalars = ['CELLS', 'material']
iPoints.ThresholdRange = [2.0, 2.0]
RenameSource('iPoints', iPoints)
# show data in view
iPointsDisplay = Show(iPoints, renderView1, 'UnstructuredGridRepresentation')
iPointsDisplay.SetRepresentationType('Point Gaussian')
iPointsDisplay.ShaderPreset = 'Plain circle'
iPointsDisplay.AmbientColor = [1.0, 0.0, 0.0]
iPointsDisplay.DiffuseColor = [1.0, 0.0, 0.0]
iPointsDisplay.GaussianRadius = 0.03

# FRACTURE OUTLINE
outline = Threshold(Input=dim1Dfull)
outline.Scalars = ['CELLS', 'material']
outline.ThresholdRange = [2.0, 2.0]
RenameSource('fOutline', outline)
# show data in view
outlineDisplay = Show(outline, renderView1, 'UnstructuredGridRepresentation')
outlineDisplay.AmbientColor = [1.0, 0.0, 0.0]
outlineDisplay.DiffuseColor = [1.0, 0.0, 0.0]
outlineDisplay.LineWidth = 1.4

# RIBS
ribs = Threshold(Input=dim1Dfull)
ribs.Scalars = ['CELLS', 'material']
ribs.ThresholdRange = [3.0, 8.0]
RenameSource('Ribs', ribs)
# show data in view
ribsDisplay = Show(ribs, renderView1, 'UnstructuredGridRepresentation')
ribsDisplay.AmbientColor = [0.0, 0.0, 0.4980392156862745]
ribsDisplay.DiffuseColor = [0.0, 0.0, 0.4980392156862745]


# SHRINK
shrink1 = Shrink(Input=elementRange)
shrink1.ShrinkFactor = 0.82

# 2D FULL
dim2Dshrink = Threshold(Input=shrink1)
dim2Dshrink.Scalars = ['CELLS', 'Dimension']
dim2Dshrink.ThresholdRange = [2.0, 2.0]
RenameSource('2D', dim2Dshrink)


# INTACT FACES
intact = Threshold(Input=dim2Dshrink)
intact.Scalars = ['CELLS', 'material']
intact.ThresholdRange = [1.0, 1.0]
RenameSource('intactFaces', intact)
# show data in view
intactDisplay = Show(intact, renderView1, 'UnstructuredGridRepresentation')
intactDisplay.AmbientColor = [0.0, 0.0, 0.4980392156862745]
intactDisplay.DiffuseColor = [0.0, 0.0, 0.4980392156862745]

# REFINEMENTS
refinements = Threshold(Input=dim2Dshrink)
refinements.Scalars = ['CELLS', 'material']
refinements.ThresholdRange = [3.0, 8.0]
RenameSource('Refinements', refinements)
# show data in view
refinementsDisplay = Show(refinements, renderView1, 'UnstructuredGridRepresentation')
# set scalar coloring
ColorBy(refinementsDisplay, ('CELLS', 'material'))
# rescale color and/or opacity maps used to include current data range
refinementsDisplay.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
refinementsDisplay.SetScalarBarVisibility(renderView1, True)
# get color transfer function/color map for 'material'
materialLUT = GetColorTransferFunction('material')
# get opacity transfer function/opacity map for 'material'
materialPWF = GetOpacityTransferFunction('material')
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
materialLUT.ApplyPreset('Blue - Green - Orange', True)
# Rescale transfer function
materialLUT.RescaleTransferFunction(3.0, 8.0)
# Rescale transfer function
materialPWF.RescaleTransferFunction(3.0, 8.0)
# hide color bar/color legend
refinementsDisplay.SetScalarBarVisibility(renderView1, False)


# SURFACE MESH
surfaceMesh = Threshold(Input=dim2Dshrink)
surfaceMesh.Scalars = ['CELLS', 'material']
surfaceMesh.ThresholdRange = [2.0, 2.0]
RenameSource('SurfaceMesh', surfaceMesh)
# show data in view
surfaceMeshDisplay = Show(surfaceMesh, renderView1, 'UnstructuredGridRepresentation')
surfaceMeshDisplay.AmbientColor = [1.0, 0.3333333333333333, 0.0]
surfaceMeshDisplay.DiffuseColor = [1.0, 0.3333333333333333, 0.0]


# 3D SHRINK
dim3Dshrink = Threshold(Input=shrink1)
dim3Dshrink.Scalars = ['CELLS', 'Dimension']
dim3Dshrink.ThresholdRange = [2.0, 3.0]
RenameSource('3D', dim3Dshrink)

# SUBMESH
submesh = Threshold(Input=dim3Dshrink)
submesh.Scalars = ['CELLS', 'material']
submesh.ThresholdRange = [3.0, 3.0]
RenameSource('SubMesh', submesh)
# show data in view
submeshDisplay = Show(submesh, renderView1, 'UnstructuredGridRepresentation')


# 2D FULL
dim2Dfull = Threshold(Input=elementRange)
dim2Dfull.Scalars = ['CELLS', 'Dimension']
dim2Dfull.ThresholdRange = [2.0, 2.0]
RenameSource('2D', dim2Dfull)

# FRACTURES FULL
Fractures = Threshold(Input=dim2Dfull)
Fractures.Scalars = ['CELLS', 'material']
Fractures.ThresholdRange = [100.0, 100.0]
RenameSource('Fractures', Fractures)
FracturesDisplay = Show(Fractures, renderView1,'UnstructuredGridRepresentation')
FracturesDisplay.SetRepresentationType('Surface With Edges')
FracturesDisplay.EdgeColor = [1.0, 1.0, 1.0]
FracturesDisplay.Opacity = 0.5

# update the view to ensure updated data information
Hide(vtkmesh, renderView1)
Hide(elementRange, renderView1)
Hide(dim3Dfull, renderView1)
Hide(dim1Dfull, renderView1)
Hide(shrink1, renderView1)
Hide(ribs, renderView1)
Hide(dim2Dshrink, renderView1)
Hide(dim3Dshrink, renderView1)
Hide(submesh, renderView1)
Hide(dim2Dfull, renderView1)
Hide(Fractures, renderView1)
renderView1.Update()


# current camera placement for renderView1
renderView1.InteractionMode = '3D'
renderView1.CameraPosition = [3.4700000286102295, 4.0, 10000.0]
renderView1.CameraFocalPoint = [3.4700000286102295, 4.0, 0.0]
renderView1.CameraParallelScale = 7.1307240904601645


# reset view to fit data
renderView1.ResetCamera()