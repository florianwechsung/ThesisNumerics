from paraview.simple import *
import os

directories = [
    "levelset-base-elasticity-cr-0.01/",
    "levelset-base-elasticity-cr-None/",
    "levelset-base-laplace-cr-0.01/",
    "levelset-base-laplace-cr-None/",
]
hist_titles = [
    r'$(u, v) = \frac{1}{\alpha} (\mathcal{B}u,\mathcal{B}v)_{L^2} + (\mathrm{sym}(\nabla u), \mathrm{sym}(\nabla v))_{L^2}$',
    r'$(u, v) = (\mathrm{sym}(\nabla u), \mathrm{sym}(\nabla v))_{L^2}$',
    r'$(u, v) = \frac{1}{\alpha} (\mathcal{B}u,\mathcal{B}v)_{L^2} + (\nabla u, \nabla v)_{L^2}$',
    r'$(u, v) = (\nabla u, \nabla v)_{L^2}$',
]
filename = "domain.pvd"

i = 0
directory = directories[i]
hist_title = hist_titles[i]
hist_maximum = 4000.
#     directory = r'/Users/florianwechsung/Documents/Uni/DPhil/Conformal-Mapping-Shape-Optimization/code/examples/levelset_fem_output/'

# filename = r'dilatation_CR+sym(grad),1e-2.pvd'
# hist_title = 

# filename = r'dilatation_CR+grad,1e-2.pvd'
# hist_title = '$(u, v) = \\frac{1}{\\alpha} (\\mathcal{B}u,\\mathcal{B}v)_{L^2} + (\\nabla u, \\nabla v)_{L^2}$'

# filename = r'dilatation_grad.pvd'
# hist_title = '$(u, v) =  (\\nabla u, \\nabla v)_{L^2}$'

# filename = r'dilatation_sym(grad).pvd'
# hist_title = '$(u, v) =  (\\mathrm{sym}(\\nabla u), \\mathrm{sym}(\\nabla v))_{L^2}$'

# filename = r'dilatation_CR+sym(grad),3e-1.pvd'
# hist_title = '$(u, v) = \\frac{1}{\\alpha} (\\mathcal{B}u,\\mathcal{B}v)_{L^2} + (\\mathrm{sym}(\\nabla u), \\mathrm{sym}(\\nabla v))_{L^2}$'

# filename = r'dilatation_CR+grad,3e-1.pvd'
# hist_title = '$(u, v) = \\frac{1}{\\alpha} (\\mathcal{B}u,\\mathcal{B}v)_{L^2} + (\\nabla u, \\nabla v)_{L^2}$'

label = filename[11:-4]

img_directory = directory + "img/"
if not os.path.exists(img_directory):
    os.makedirs(img_directory)
# create a new 'PVD Reader'
dilatation_pvd = PVDReader(FileName=directory + filename)
# get animation scene
animationScene = GetAnimationScene()

# update animation scene based on data timesteps
animationScene.UpdateAnimationUsingDataTimeSteps()
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.OrientationAxesVisibility = 0
renderView1.Background = [1, 1, 1]
dilatation_display = Show(dilatation_pvd, renderView1)
dilatation_display.Representation = 'Surface'
dilatation_display.ColorArrayName = [None, '']
# dilatation_display.OSPRayScaleArray = 'function_16'
# dilatation_display.OSPRayScaleFunction = 'PiecewiseFunction'
dilatation_display.SelectOrientationVectors = 'None'
dilatation_display.SelectScaleArray = 'None'
dilatation_display.GlyphType = 'Arrow'
# dilatation_display.PolarAxes = 'PolarAxesRepresentation'
dilatation_display.ScalarOpacityUnitDistance = 0.417380156334143
# dilatation_display.GaussianRadius = 0.5
# dilatation_display.SetScaleArray = ['POINTS', 'function_16']
# dilatation_display.ScaleTransferFunction = 'PiecewiseFunction'
# dilatation_display.OpacityArray = ['POINTS', 'function_16']
# dilatation_display.OpacityTransferFunction = 'PiecewiseFunction'

# change representation type
dilatation_display.SetRepresentationType('Surface With Edges')

animationScene.GoToLast()
renderView1.ViewSize = [500, 500]

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [-0.7, 0.8, 1]
renderView1.CameraFocalPoint = [-0.7, 0.8, 0]
renderView1.CameraParallelScale = 1.15
# SaveScreenshot(directory + label + "_mesh.png", renderView1, ImageResolution=(1000, 1000))
SaveScreenshot(img_directory + label + "_mesh.png", view=renderView1, magnification=2)
renderView1.CameraPosition = [-0.73, 0.73, 1]
renderView1.CameraFocalPoint = [-0.73, 0.73, 0]
renderView1.CameraParallelScale = 0.06
# SaveScreenshot(directory + label + "_mesh_zoom.png", renderView1, ImageResolution=(1000, 1000))
SaveScreenshot(img_directory + label + "_mesh_zoom.png", view=renderView1, magnification=2)

renderView1.ResetCamera()

# #changing interaction mode based on data extents
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [0.0, 0.0, 10000.0]

# # update the view to ensure updated data information
renderView1.Update()

# # change representation type
# dilatation_pvd.SetRepresentationType('Surface With Edges')

# # Hide orientation axes
# renderView1.OrientationAxesVisibility = 0

# create a new 'Mesh Quality'
meshQuality1 = MeshQuality(Input=dilatation_pvd)

# get color transfer function/color map for 'Quality'
qualityLUT = GetColorTransferFunction('Quality')

# get opacity transfer function/opacity map for 'Quality'
qualityPWF = GetOpacityTransferFunction('Quality')

# show data in view
meshQuality1Display = Show(meshQuality1, renderView1)
# trace defaults for the display properties.
meshQuality1Display.Representation = 'Surface'
meshQuality1Display.AmbientColor = [0.0, 0.0, 0.0]
meshQuality1Display.ColorArrayName = ['CELLS', 'Quality']
meshQuality1Display.LookupTable = qualityLUT
# meshQuality1Display.OSPRayScaleArray = 'Quality'
# meshQuality1Display.OSPRayScaleFunction = 'PiecewiseFunction'
meshQuality1Display.SelectOrientationVectors = 'None'
meshQuality1Display.ScaleFactor = 0.6000000000000001
meshQuality1Display.SelectScaleArray = 'Quality'
meshQuality1Display.GlyphType = 'Arrow'
# meshQuality1Display.GlyphTableIndexArray = 'Quality'
# meshQuality1Display.DataAxesGrid = 'GridAxesRepresentation'
# meshQuality1Display.PolarAxes = 'PolarAxesRepresentation'
meshQuality1Display.ScalarOpacityFunction = qualityPWF
meshQuality1Display.ScalarOpacityUnitDistance = 0.30725553282259715
# meshQuality1Display.GaussianRadius = 0.30000000000000004
# meshQuality1Display.SetScaleArray = ['POINTS', 'function_10']
# meshQuality1Display.ScaleTransferFunction = 'PiecewiseFunction'
# meshQuality1Display.OpacityArray = ['POINTS', 'function_10']
# meshQuality1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
# meshQuality1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
# meshQuality1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
# meshQuality1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
# meshQuality1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
# meshQuality1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
# meshQuality1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# hide data in view
Hide(dilatation_pvd, renderView1)

# show color bar/color legend
meshQuality1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# CreateLayout('Layout #2')


# set active view
SetActiveView(None)

# Create a new 'Histogram View'
histogramView1 = CreateView('XYHistogramChartView')
histogramView1.ViewSize = [400, 200]
histogramView1.LeftAxisRangeMaximum = 6.66
histogramView1.BottomAxisRangeMaximum = 6.66
histogramView1.ChartTitle = hist_title
histogramView1.BottomAxisTitle = "Mesh quality"
histogramView1.LeftAxisUseCustomRange = 1
histogramView1.LeftAxisRangeMaximum = hist_maximum
histogramView1.LeftAxisRangeMinimum = 0.0
# get layout
layout2 = GetLayout()

# place view in the layout
layout2.AssignView(0, histogramView1)

# set active source
SetActiveSource(meshQuality1)

# show data in view
meshQuality1Display_1 = Show(meshQuality1, histogramView1)
# trace defaults for the display properties.
meshQuality1Display_1.SelectInputArray = ['CELLS', 'Quality']
meshQuality1Display_1.UseCustomBinRanges = 1
meshQuality1Display_1.CustomBinRanges = [1+1e-5, 2-1e-5]

# Properties modified on meshQuality1Display_1
meshQuality1Display_1.BinCount = 40

### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 1.0

# animationScene.GoToLast()
ExportView(img_directory + label + "_mesh_quality_final.eps", view=histogramView1)
animationScene.GoToFirst()
histogramView1.ChartTitle = "Initial mesh quality"
ExportView(img_directory + label + "_mesh_quality_initial.eps", view=histogramView1)
