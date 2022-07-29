# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

from nc_reader import nc_reader
import numpy as np
import os
import colorcet as cc

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

class iso_surface:

    def __init__(self):
        self._ncreader = nc_reader()
        self._pvnc = None

        # find settings proxy
        colorPalette = GetSettingsProxy('ColorPalette')
        self.colormap = 'Cool to Warm'

        self._field_label = {
            'vorticity_magnitude': 'vorticity magnitude',
            'helicity': 'helicity'
        }


        # black font
        colorPalette.Text = [0.0, 0.0, 0.0]
        colorPalette.Foreground = [0.0, 0.0, 0.0]

        create_cmap = False
        if create_cmap:
            self._create_colormap()

    def _create_colormap(self):
        """
        27 July 2022
        https://discourse.paraview.org/t/how-to-export-paraview-colormap-into-a-format-that-could-be-read-by-matplotlib/2436/2
        https://www.paraview.org/Wiki/Colormaps
        """
        cmap = cc.cm['rainbow4']
        scheme = 'rainbow4'

        with open('rainbow4.xml', 'w') as fid:
            fid.write('<ColorMaps>\n')
            fid.write('<ColorMap name="{}" space="RGB">\n'.format(scheme))
            N = cmap.N
            scale = np.linspace(0, 1, cmap.N)
            colors = cmap(scale)[:, :]
            for i in range(N):
                x = [(i-1)/(N-1)] + colors[i,:].tolist()
                fid.write(
                    '<Point x="{:2f}" o="{:2f}" r="{:2f}" g="{:2f}" b="{:2f}" />\n'.format(x[0],x[4],
                                                                                           x[1],
                                                                                           x[2],
                                                                                           x[3]))

            fid.write('</ColorMap>\n')
            fid.write('</ColorMaps>')
        ImportPresets(filename='rainbow4.xml')

    def open(self, fname, **kwargs):
        self._ncreader.open(fname)
        basename = os.path.basename(fname)
        self._times = self._ncreader.get_all('t')
        self._nsteps = len(self._times)
        self._pvnc = NetCDFReader(registrationName=basename, FileName=[fname])
        self._pvnc.Dimensions = '(z, y, x)'
        self._pvnc.SphericalCoordinates = 0
        self._animation_scene = GetAnimationScene()
        self._animation_scene.UpdateAnimationUsingDataTimeSteps()
        self._set_basic_render_view()
        self._create_time_stamp_filter()
        self._create_programmable_filters()

        self._layout = GetLayout()
        # layout/tab size in pixels
        self._width = kwargs.pop('width', 1951)
        self._height = kwargs.pop('height', 1660)
        self._layout.SetSize(self._width, self._height)

    def render(self, field_name, step, n_iso, **kwargs):
        field_data = self._ncreader.get_dataset(step=step, name=field_name)
        field_data = field_data ** 2
        vmax = kwargs.pop('vmax', field_data.max())
        vmin = kwargs.pop('vmin', field_data.min())

        self._animation_scene.AnimationTime = self._times[step]
        self._create_contours(field_name, vmin=vmin, vmax=vmax, n_iso=n_iso)
        self._create_color_bar(field_name=field_name, vmax=vmax)
        self._create_surface(field_name)
        self._set_camera_position()

    def save_camera_orbiting_animation(self, field_name, step, n_frames, **kwargs):
        """
        29 July 2022
        https://discourse.paraview.org/t/animation-camera-orbit-python/2907/3
        """
        tmp_dir = kwargs.pop('tmp_dir', 'temp_dir')
        n_iso = kwargs.pop('n_iso', 20)

        if os.path.exists(tmp_dir):
            print("Error: Directory '" + tmp_dir + "' already exists. Exiting.")
            exit()

        os.mkdir(tmp_dir)

        file_name = kwargs.pop('file_name', 'orbit_movie.mp4')
        file_path = kwargs.pop('file_path', './')
        fps = kwargs.pop('fps', 25)
        keep_frames = kwargs.pop('keep_frames', False)

        self.render(field_name=field_name, step=step, n_iso=n_iso)

        self._render_view.CameraPosition = [-10, 4, 4]
        self._render_view.CameraViewUp = [0.0, 0.0, 1.0]
        self._render_view.CameraFocalPoint = [0, 0, 0]
        self._render_view.AxesGrid.Visibility = 0

        self._color_bar.ScalarBarLength = 0.25
        self._color_bar.Position = [0.89, 0.15]

        self._render_view.Update()

        camera = self._render_view.GetActiveCamera()

        dtheta = 360 / n_frames
        for i in range(0, n_frames):
            camera.Azimuth(dtheta)
            self._render_view.Update()
            self.export(file_path=tmp_dir, file_name='frame' + str(i).zfill(5) + '.png')

        os.system('ffmpeg -i ' + os.path.join(tmp_dir, 'frame%05d.png') +
                  ' -c:v libx264 -vf fps=' + str(fps) + ' ' + file_name)

        if not keep_frames:
            for i in range(0, n_frames):
                os.remove(os.path.join(tmp_dir, 'frame' + str(i).zfill(5) + '.png'))
            os.rmdir(tmp_dir)

    def export(self, file_path, file_name):
        # make sure we have recent view
        self._render_view.Update()
        # get extension
        _, ext = os.path.splitext(file_name)
        if ext == '.png':
            SaveScreenshot(os.path.join(file_path, file_name),
                           self._render_view,
                           ImageResolution=[self._width, self._height],
                           CompressionLevel=5)
        elif ext == '.jpg':
            SaveScreenshot(os.path.join(file_path, file_name),
                           self._render_view,
                           ImageResolution=[self._width, self._height],
                           Quality=50)
        elif ext == '.eps':
            ExportView(os.path.join(file_path, file_name), view=self._render_view)
        else:
            print("The " + ext.upper() + " file format is not supported.")

    def save_animation(self, field_name, beg, end, **kwargs):

        tmp_dir = kwargs.pop('tmp_dir', 'temp_dir')
        n_iso = kwargs.pop('n_iso', 20)

        if os.path.exists(tmp_dir):
            print("Error: Directory '" + tmp_dir + "' already exists. Exiting.")
            exit()

        os.mkdir(tmp_dir)

        file_name = kwargs.pop('file_name', 'beltrami_instability.mp4')
        file_path = kwargs.pop('file_path', './')
        fps = kwargs.pop('fps', 25)
        keep_frames = kwargs.pop('keep_frames', False)


        for i in range(beg, end+1):
            self.render(field_name=field_name, step=i, n_iso=n_iso, vmin=0.0)
            self.export(file_path=tmp_dir, file_name='frame' + str(i).zfill(5) + '.png')
            self._clear()

        #SaveAnimation(os.path.join(tmp_dir, 'frame.png'),
                      #self._render_view, ImageResolution=[self._width, self._height],
                      #FrameWindow=[beg, end])

        os.system('ffmpeg -r 2 -i ' + os.path.join(tmp_dir, 'frame%05d.png') +
                  ' -c:v libx264 -vf fps=' + str(fps) + ' ' + file_name)

        if not keep_frames:
            for i in range(beg, end+1):
                os.remove(os.path.join(tmp_dir, 'frame' + str(i).zfill(5) + '.png'))
            os.rmdir(tmp_dir)

    def _clear(self):
        self._contour_display.Visibility = 0
        self._color_bar.Visibility = 0
        self._prog_filter_display.Visibility = 0
        self._render_view.Update()

    def close(self):
        self._ncreader.close()
        self._pvnc = None

    def _set_basic_render_view(self):
        self._render_view = GetActiveViewOrCreate('RenderView')
        self._render_view.UseColorPaletteForBackground = 0
        self._render_view.Background = [1, 1, 1] # white background
        self._render_view.ResetCamera(False)
        self._render_view.CenterAxesVisibility = 0
        self._render_view.OrientationAxesVisibility = 0
        self._render_view.AxesGrid.Visibility = 1
        self._render_view.AxesGrid.DataBoundsScaleFactor = 1.008

        axis_labels = [-1.5, -0.75, 0.0, 0.75, 1.5]
        self._render_view.AxesGrid.XAxisUseCustomLabels = 1
        self._render_view.AxesGrid.XAxisLabels = axis_labels[0:-1]

        self._render_view.AxesGrid.YAxisUseCustomLabels = 1
        self._render_view.AxesGrid.YAxisLabels = axis_labels

        self._render_view.AxesGrid.ZAxisUseCustomLabels = 1
        self._render_view.AxesGrid.ZAxisLabels = axis_labels

        self._render_view.AxesGrid.XTitle = 'x'
        self._render_view.AxesGrid.YTitle = 'y'
        self._render_view.AxesGrid.ZTitle = 'z'
        self._render_view.AxesGrid.XTitleFontFamily = 'Courier'
        self._render_view.AxesGrid.XTitleFontSize = 30
        self._render_view.AxesGrid.YTitleFontFamily = 'Courier'
        self._render_view.AxesGrid.YTitleFontSize = 30
        self._render_view.AxesGrid.ZTitleFontFamily = 'Courier'
        self._render_view.AxesGrid.ZTitleFontSize = 30
        self._render_view.AxesGrid.XLabelFontFamily = 'Courier'
        self._render_view.AxesGrid.XLabelFontSize = 25
        self._render_view.AxesGrid.YLabelFontFamily = 'Courier'
        self._render_view.AxesGrid.YLabelFontSize = 25
        self._render_view.AxesGrid.ZLabelFontFamily = 'Courier'
        self._render_view.AxesGrid.ZLabelFontSize = 25
        self._render_view.Update()

    def _create_time_stamp_filter(self):
        time_filter = AnnotateTimeFilter(registrationName='TimeStampFilter', Input=self._pvnc)
        time_filter.Format = 't = {time:3.5f}'

        SetActiveSource(time_filter)

        time_filter_display = Show(time_filter, self._render_view, 'TextSourceRepresentation')
        time_filter_display.FontFamily = 'Courier'
        time_filter_display.FontSize = 40
        time_filter_display.Bold = 1
        time_filter_display.WindowLocation = 'Any Location'
        time_filter_display.Position = [0.05, 0.9]
        self._render_view.Update()

    def _create_programmable_filters(self):
        #
        # vorticity magnitude
        #
        self._prog_filter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=self._pvnc)
        self._prog_filter1.Script = ''
        self._prog_filter1.RequestInformationScript = ''
        self._prog_filter1.RequestUpdateExtentScript = ''
        self._prog_filter1.PythonPath = ''

        self._prog_filter1.Script = """
import numpy as np
xi = inputs[0].PointData['x_vorticity']
eta = inputs[0].PointData['y_vorticity']
zeta = inputs[0].PointData['z_vorticity']
output.PointData.append(np.sqrt(xi ** 2 + eta ** 2 + zeta ** 2), 'vorticity_magnitude')"""
        self._prog_filter1.RequestInformationScript = ''
        self._prog_filter1.RequestUpdateExtentScript = ''
        self._prog_filter1.PythonPath = ''

        #
        # helicity
        #
        self._prog_filter2 = ProgrammableFilter(registrationName='ProgrammableFilter2', Input=self._pvnc)
        self._prog_filter2.Script = ''
        self._prog_filter2.RequestInformationScript = ''
        self._prog_filter2.RequestUpdateExtentScript = ''
        self._prog_filter2.PythonPath = ''

        self._prog_filter2.Script = """
import numpy as np
u = inputs[0].PointData['x_velocity']
v = inputs[0].PointData['y_velocity']
w = inputs[0].PointData['z_velocity']
xi = inputs[0].PointData['x_vorticity']
eta = inputs[0].PointData['y_vorticity']
zeta = inputs[0].PointData['z_vorticity']
output.PointData.append(u * xi + v * eta + w * zeta, 'helicity')"""
        self._prog_filter2.RequestInformationScript = ''
        self._prog_filter2.RequestUpdateExtentScript = ''
        self._prog_filter2.PythonPath = ''

        self._prog_filters = {
            'vorticity_magnitude': self._prog_filter1,
            'helicity': self._prog_filter2
        }

    def _create_contours(self, field_name, vmin, vmax, n_iso):
        contour = Contour(registrationName='Contour1', Input=self._prog_filters[field_name])
        contour.ContourBy = ['POINTS', field_name]
        contour.Isosurfaces = np.linspace(vmin, vmax, n_iso)
        contour.PointMergeMethod = 'Uniform Binning'

        # set active source
        SetActiveSource(contour)

        # show data in view
        self._contour_display = Show(contour, self._render_view, 'GeometryRepresentation')

        # set scalar coloring
        ColorBy(self._contour_display, ('POINTS', field_name))

        # rescale color and/or opacity maps used to include current data range
        self._contour_display.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        self._contour_display.SetScalarBarVisibility(self._render_view, True)

        # get color transfer function/color map for field_name
        self._lut = GetColorTransferFunction(field_name)

        self._lut.EnableOpacityMapping = 1

        self._lut.ApplyPreset(self.colormap, True)

        # get opacity transfer function/opacity map for field_name
        self._pwf = GetOpacityTransferFunction(field_name)

        ## Rescale transfer function
        self._lut.RescaleTransferFunction(0.0, vmax)

        ## Rescale transfer function
        self._pwf.RescaleTransferFunction(0.0, vmax)

        # trace defaults for the display properties.
        self._contour_display.Representation = 'Surface'
        self._contour_display.ColorArrayName = ['POINTS', field_name]
        self._contour_display.SelectTCoordArray = 'None'
        self._contour_display.SelectNormalArray = 'Normals'
        self._contour_display.SelectTangentArray = 'None'
        self._contour_display.OSPRayScaleArray = field_name
        self._contour_display.OSPRayScaleFunction = 'PiecewiseFunction'
        self._contour_display.SelectOrientationVectors = 'None'
        self._contour_display.ScaleFactor = 0.30434179306030273
        self._contour_display.SelectScaleArray = field_name
        self._contour_display.GlyphType = 'Arrow'
        self._contour_display.GlyphTableIndexArray = field_name
        self._contour_display.GaussianRadius = 0.015217089653015136
        self._contour_display.SetScaleArray = ['POINTS', field_name]
        self._contour_display.ScaleTransferFunction = 'PiecewiseFunction'
        self._contour_display.OpacityArray = ['POINTS', field_name]
        self._contour_display.OpacityTransferFunction = 'PiecewiseFunction'
        self._contour_display.DataAxesGrid = 'GridAxesRepresentation'
        self._contour_display.PolarAxes = 'PolarAxesRepresentation'
        self._contour_display.LookupTable = self._lut
        self._render_view.Update()

    def _create_color_bar(self, field_name, vmax):
        self._color_bar = GetScalarBar(self._lut, self._render_view)

        self._color_bar.TitleFontFamily = 'Courier'
        self._color_bar.TitleFontSize = 30
        self._color_bar.LabelFontFamily = 'Courier'
        self._color_bar.LabelFontSize = 30

        self._color_bar.ScalarBarThickness = 20
        self._color_bar.ScalarBarLength = 0.5
        self._color_bar.WindowLocation = 'Any Location'
        self._color_bar.Position = [0.85, 0.25]
        self._color_bar.Title = self._field_label[field_name]

        self._color_bar.TitleJustification = 'Centered'

        self._color_bar.AddRangeLabels = 0
        self._color_bar.RangeLabelFormat = ''

        self._color_bar.AutomaticLabelFormat = 0
        self._color_bar.LabelFormat = '%-#6.2g'

        self._color_bar.UseCustomLabels = 1
        self._color_bar.CustomLabels = np.linspace(0, vmax, 10)
        self._render_view.Update()

    def _create_surface(self, field_name):
        self._prog_filter_display = Show(self._prog_filters[field_name],
                                         self._render_view, 'UniformGridRepresentation')

        SetActiveSource(self._prog_filter_display)

        self._prog_filter_display.SetRepresentationType('Surface')
        self._prog_filter_display.Opacity = 0.5
        self._prog_filter_display.ScaleFactor = 0.323976504603413
        self._prog_filter_display.DataAxesGrid = 'GridAxesRepresentation'
        self._prog_filter_display.PolarAxes = 'PolarAxesRepresentation'
        self._prog_filter_display.SetScaleArray = ['POINTS', field_name]
        self._prog_filter_display.ScaleTransferFunction = 'PiecewiseFunction'
        self._prog_filter_display.ColorArrayName = ['POINTS', field_name]
        self._prog_filter_display.OpacityArray = ['POINTS', field_name]
        self._prog_filter_display.OpacityTransferFunction = 'PiecewiseFunction'
        self._prog_filter_display.LookupTable = self._lut
        self._render_view.Update()

    def _set_camera_position(self):
        self._render_view.CameraPosition = [-10, 4, 4]
        self._render_view.CameraViewUp = [0.0, 0.0, 1.0]
        self._render_view.CameraFocalPoint = [-0.6, 0, 0]
        self._render_view.CameraParallelScale = 2.0
        self._render_view.Update()
