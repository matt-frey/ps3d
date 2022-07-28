# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

from nc_reader import nc_reader
import numpy as np

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find settings proxy
colorPalette = GetSettingsProxy('ColorPalette')

# Properties modified on colorPalette --> black font
colorPalette.Text = [0.0, 0.0, 0.0]

fname = '/home/matthias/Documents/projects/ps3d/examples/beltrami_32_fields.nc'
step = 60
niso = 10

ncreader = nc_reader()
ncreader.open(fname)
times = ncreader.get_all('t')
vmag = ncreader.get_dataset(step=step, name='vorticity_magnitude')
vmin = vmag.min()
vmax = vmag.max()
ncreader.close()

print(vmin, vmax)

# create a new 'NetCDF Reader'
beltrami_32_fieldsnc = NetCDFReader(
    registrationName='beltrami_32_fields.nc',
    FileName=[fname])

beltrami_32_fieldsnc.Dimensions = '(z, y, x)'

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on animationScene1
animationScene1.AnimationTime = times[step]

# Properties modified on beltrami_32_fieldsnc
beltrami_32_fieldsnc.SphericalCoordinates = 0

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.UseColorPaletteForBackground = 0

# Properties modified on renderView1
renderView1.Background = [1, 1, 1]

renderView1.ResetCamera(False)

renderView1.CenterAxesVisibility = 1

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.DataBoundsScaleFactor = 1.008

axis_labels = [-1.5, -0.75, 0.0, 0.75, 1.5]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = axis_labels[0:-1]

renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = axis_labels

renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = axis_labels

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'x'
renderView1.AxesGrid.YTitle = 'y'
renderView1.AxesGrid.ZTitle = 'z'
renderView1.AxesGrid.XTitleFontFamily = 'Courier'
renderView1.AxesGrid.XTitleFontSize = 30
renderView1.AxesGrid.YTitleFontFamily = 'Courier'
renderView1.AxesGrid.YTitleFontSize = 30
renderView1.AxesGrid.ZTitleFontFamily = 'Courier'
renderView1.AxesGrid.ZTitleFontSize = 30
renderView1.AxesGrid.XLabelFontFamily = 'Courier'
renderView1.AxesGrid.XLabelFontSize = 25
renderView1.AxesGrid.YLabelFontFamily = 'Courier'
renderView1.AxesGrid.YLabelFontSize = 25
renderView1.AxesGrid.ZLabelFontFamily = 'Courier'
renderView1.AxesGrid.ZLabelFontSize = 25

renderView1.Update()

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=beltrami_32_fieldsnc)
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
programmableFilter1.Script = """import numpy as np
xvor = inputs[0].PointData['x_vorticity']
yvor = inputs[0].PointData['y_vorticity']
zvor = inputs[0].PointData['z_vorticity']
output.PointData.append(xvor ** 2 + yvor ** 2 + zvor ** 2, 'vorticity_magnitude')"""
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''


############################################################################################################
# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=programmableFilter1)
contour1.ContourBy = ['POINTS', 'vorticity_magnitude']
#contour1.Isosurfaces = [4.497020606762032]
contour1.Isosurfaces = np.linspace(vmin, vmax, niso)
contour1.PointMergeMethod = 'Uniform Binning'

# set active source
SetActiveSource(contour1)

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# set scalar coloring
ColorBy(contour1Display, ('POINTS', 'vorticity_magnitude'))

#https://docs.paraview.org/en/latest/ReferenceManual/colorMapping.html

# rescale color and/or opacity maps used to include current data range
contour1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vorticity_magnitude'
vorticity_magnitudeLUT = GetColorTransferFunction('vorticity_magnitude')

vorticity_magnitudeLUT.EnableOpacityMapping = 1

vorticity_magnitudeLUT.ApplyPreset('MyRainbow', True)
#vorticity_magnitudeLUT.ApplyPreset('RedOpaqueBlue', True)

# get opacity transfer function/opacity map for 'vorticity_magnitude'
vorticity_magnitudePWF = GetOpacityTransferFunction('vorticity_magnitude')

## Rescale transfer function
vorticity_magnitudeLUT.RescaleTransferFunction(0.0, vmax)

## Rescale transfer function
vorticity_magnitudePWF.RescaleTransferFunction(0.0, vmax)

# Properties modified on vorticity_magnitudePWF
#vorticity_magnitudePWF.Points = np.linspace(0, 1, niso) #[vmin, vmax]

# Properties modified on vorticity_magnitudeLUT
#vorticity_magnitudeLUT.RGBPoints = np.linspace(0, 1, niso)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'vorticity_magnitude']
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'vorticity_magnitude'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.30434179306030273
contour1Display.SelectScaleArray = 'vorticity_magnitude'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'vorticity_magnitude'
contour1Display.GaussianRadius = 0.015217089653015136
contour1Display.SetScaleArray = ['POINTS', 'vorticity_magnitude']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'vorticity_magnitude']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'
contour1Display.LookupTable = vorticity_magnitudeLUT

# get color legend/bar for vorticity_magnitudeLUT in view renderView1
vmagLUTColorBar = GetScalarBar(vorticity_magnitudeLUT, renderView1)

# Properties modified on pressureLUTColorBar
vmagLUTColorBar.TitleFontFamily = 'Courier'
vmagLUTColorBar.TitleFontSize = 30
vmagLUTColorBar.LabelFontFamily = 'Courier'
vmagLUTColorBar.LabelFontSize = 30

# Properties modified on pressureLUTColorBar
vmagLUTColorBar.ScalarBarThickness = 20
vmagLUTColorBar.ScalarBarLength = 0.5
vmagLUTColorBar.WindowLocation = 'Any Location'
vmagLUTColorBar.Position = [0.85, 0.25]
vmagLUTColorBar.Title = 'vorticity magnitude squared'

# Properties modified on pressureLUTColorBar
vmagLUTColorBar.TitleJustification = 'Centered'


vmagLUTColorBar.AddRangeLabels = 0
vmagLUTColorBar.RangeLabelFormat = ''

vmagLUTColorBar.AutomaticLabelFormat = 0
vmagLUTColorBar.LabelFormat = '%-#6.2g'

# Properties modified on vorticity_magnitudeLUTColorBar
vmagLUTColorBar.UseCustomLabels = 1
vmagLUTColorBar.CustomLabels = np.linspace(0, vmax, 10)


renderView1.Update()

###########################################################################################

# show data in view
programmableFilter1Display = Show(programmableFilter1, renderView1, 'UniformGridRepresentation')

SetActiveSource(programmableFilter1Display)

# surface:
#beltrami_32_fieldsncDisplay = Show(beltrami_32_fieldsnc, renderView1, 'StructuredGridRepresentation')
programmableFilter1Display.SetRepresentationType('Surface')
programmableFilter1Display.Opacity = 0.5
programmableFilter1Display.ScaleFactor = 0.323976504603413
programmableFilter1Display.DataAxesGrid = 'GridAxesRepresentation'
programmableFilter1Display.PolarAxes = 'PolarAxesRepresentation'
programmableFilter1Display.SetScaleArray = ['POINTS', 'vorticity_magnitude']
programmableFilter1Display.ScaleTransferFunction = 'PiecewiseFunction'
programmableFilter1Display.ColorArrayName = ['POINTS', 'vorticity_magnitude']
programmableFilter1Display.OpacityArray = ['POINTS', 'vorticity_magnitude']
programmableFilter1Display.OpacityTransferFunction = 'PiecewiseFunction'
programmableFilter1Display.LookupTable = vorticity_magnitudeLUT


## Rescale transfer function
#vorticity_magnitudeLUT.RescaleTransferFunction(1.870760725235395e-16, 8.99365512076521)

## Rescale transfer function
#vorticity_magnitudePWF.RescaleTransferFunction(1.870760725235395e-16, 8.99365512076521)

# Properties modified on vorticity_magnitudePWF
#vorticity_magnitudePWF.Points = [1.870760725235395e-16, 1.0, 0.5, 0.0, 1.5349522829055786, 0.0, 0.5, 0.0, 3.181980515339463, 1.0, 0.5, 0.0]

# Properties modified on vorticity_magnitudeLUT
#vorticity_magnitudeLUT.RGBPoints = [1.870760725235395e-16, 0.231373, 0.298039, 0.752941, 1.5739352703094482, 0.865003, 0.865003, 0.865003, 3.181980515339463, 0.705882, 0.0156863, 0.14902]


renderView1.Update()




# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1951, 1660)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
#renderView1.CameraPosition = [-9.630692177788971, 3.56240002835346, 1.0579772867006316]
#renderView1.CameraViewUp = [0.011834494288939401, -0.2520690992266268, 0.9676368709180123]
renderView1.CameraPosition = [-10, 4, 4]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalPoint = [-0.6, 0, 0] #[-0.04908740520477295, -0.04908740520477295, -1.1192220824532735e-18]
renderView1.CameraParallelScale = 2.0 #.664319347929923

renderView1.Update()

# export view
#ExportView('/home/matthias/Documents/projects/ps3d/plotting/test_paraview.eps', view=renderView1)

# save screenshot
SaveScreenshot('/home/matthias/Documents/projects/ps3d/plotting/test_paraview.png', renderView1,
               ImageResolution=[1951, 1660])
