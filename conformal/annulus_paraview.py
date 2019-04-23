from paraview.simple import *
import os
import argparse

base = r'./output/annulus/'

parser = argparse.ArgumentParser()
parser.add_argument("--base-inner", type=str, default="elasticity",
                    choices=["elasticity", "laplace"])
parser.add_argument("--alpha", type=float, default=None)
parser.add_argument("--weighted", default=False, action="store_true")
parser.add_argument("--rstar", type=float, default=0.79)

args = parser.parse_args()
case = r"base-%s-cr-%s-weighted-%s-rstar-%.2f" % (args.base_inner, args.alpha, args.weighted, args.rstar)

# filename = r'dilatation_CR+sym(grad),1e-2,wmu,07.pvd'
filename = r'/domain.pvd'


img_directory = base + "img/"
try:
    os.makedirs(img_directory)
except:
    pass
# create a new 'PVD Reader'
dilatation_pvd = PVDReader(FileName=base + case + filename)
case = case.replace(".", "p")
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
renderView1.CameraPosition = [-0.5, 0.5, 1]
renderView1.CameraFocalPoint = [-0.5, 0.5, 0]
renderView1.CameraParallelScale = 0.5
SaveScreenshot(img_directory + case + "_mesh.png", view=renderView1, magnification=2, quality=95)
animationScene.GoToFirst()
SaveScreenshot(img_directory + case + "_mesh_initial.png", view=renderView1, magnification=2, quality=95)
