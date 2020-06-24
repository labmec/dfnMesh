
# vtkmesh
# ElementRange
#     3D
#         CoarseMesh
#     1D
#         iPoints
#         Outline
#         Ribs intact
#         Ribs cut
#     shrink
#         2D
#             intactFaces
#             Refinements
#             SurfaceMesh
#         3D
#             submesh
#     2D
#         Polygons





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

# create a new 'Threshold'
threshold1 = Threshold(Input=vtkmesh)
threshold1.Scalars = ['CELLS', 'elIndex']
threshold1.ThresholdRange = [0.0, 10000.0]
RenameSource('ElementRange', threshold1)

# create a new 'Threshold'
threshold1_1 = Threshold(Input=threshold1)
threshold1_1.Scalars = ['CELLS', 'Dimension']
threshold1_1.ThresholdRange = [3.0, 3.0]
RenameSource('3D', threshold1_1)
# create a new 'Threshold'
threshold1_1_1 = Threshold(Input=threshold1_1)
threshold1_1_1.Scalars = ['CELLS', 'material']
threshold1_1_1.ThresholdRange = [1.0, 1.0]
RenameSource('CoarseMesh', threshold1_1_1)


# create a new 'Threshold'
threshold1_2 = Threshold(Input=threshold1)
threshold1_2.Scalars = ['CELLS', 'Dimension']
threshold1_2.ThresholdRange = [1.0, 1.0]
RenameSource('1D', threshold1_2)
# create a new 'Threshold'
threshold1_2_1 = Threshold(Input=threshold1_2)
threshold1_2_1.Scalars = ['CELLS', 'material']
threshold1_2_1.ThresholdRange = [2.0, 2.0]
RenameSource('iPoints', threshold1_2_1)



#### uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).