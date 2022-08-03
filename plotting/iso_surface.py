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

    def __init__(self, create_cmaps=False):
        self._ncreader = nc_reader()
        self._pvnc = None

        # find settings proxy
        colorPalette = GetSettingsProxy('ColorPalette')

        self._derived_fields = [
            'vorticity_magnitude',
            'helicity'
        ]

        self._field_label = {
            'vorticity_magnitude': 'vorticity magnitude'
        }


        # black font
        colorPalette.Text = [0.0, 0.0, 0.0]
        colorPalette.Foreground = [0.0, 0.0, 0.0]

        if create_cmaps:
            self._create_colormaps()

    def _create_colormaps(self):
        """
        27 July 2022
        https://discourse.paraview.org/t/how-to-export-paraview-colormap-into-a-format-that-could-be-read-by-matplotlib/2436/2
        https://www.paraview.org/Wiki/Colormaps
        """
        cc_maps = ['rainbow4', 'blues', 'coolwarm']
        cc_names = ['rainbow4', 'cc_blues', 'cc_coolwarm']

        for j, name in enumerate(cc_maps):
            cmap = cc.cm[name]

            with open(cc_names[j] + '.xml', 'w') as fid:
                fid.write('<ColorMaps>\n')
                fid.write('<ColorMap name="{}" space="RGB">\n'.format(cc_names[j]))
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
                
            ImportPresets(filename=cc_names[j] + '.xml')

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

        if kwargs.get('add_time', True):
            self._create_time_stamp_filter()
        self._create_programmable_filters()

        self._layout = GetLayout()
        # layout/tab size in pixels
        self._width = kwargs.get('width', 1951)
        self._height = kwargs.get('height', 1660)
        self._layout.SetSize(self._width, self._height)

    def render(self, field_name, step, **kwargs):
        field_data = self._ncreader.get_dataset(step=step, name=field_name)
        vmax = kwargs.get('vmax', field_data.max())
        vmin = kwargs.get('vmin', field_data.min())
        if vmin is None:
            vmin = field_data.min()
        if vmax is None:
            vmax = field_data.max()        
        n_iso = kwargs.get('n_iso', 40)
        self.colormap = kwargs.get('colormap', 'Cool to Warm')
        self._invert_colormap = kwargs.get('invert_colormap', False)
        self._enable_opacity = kwargs.get('enable_opacity', False)
        self._opacity_vmin = kwargs.get('opacity_vmin', 1.0)
        self._opacity_vmax = kwargs.get('opacity_vmax', 1.0)
        self._opacity_points = kwargs.get('opacity_points', [])
        self._opacity_values = kwargs.get('opacity_values', [])
        self._color_vmin = kwargs.get('color_vmin', None)
        self._color_vmax = kwargs.get('color_vmax', None)
        self._color_points = kwargs.get('color_points', [])
        self._color_values = kwargs.get('color_values', [])
        self._add_clabel = kwargs.get('add_clabel', True)
        self._use_log_scale = kwargs.get('use_log_scale', False)
        self._n_color_bar_ticks = kwargs.get('n_color_bar_ticks', 10)

        self._animation_scene.AnimationTime = self._times[step]
        self._create_contours(field_name, vmin=vmin, vmax=vmax, n_iso=n_iso)
        self._create_color_bar(field_name=field_name, vmin=vmin, vmax=vmax)
        self._set_camera_position()

    def save_camera_orbiting_animation(self, field_name, step, n_frames, **kwargs):
        """
        29 July 2022
        https://discourse.paraview.org/t/animation-camera-orbit-python/2907/3
        """
        tmp_dir = kwargs.get('tmp_dir', 'temp_dir')

        if os.path.exists(tmp_dir):
            print("Error: Directory '" + tmp_dir + "' already exists. Exiting.")
            exit()

        os.mkdir(tmp_dir)

        file_name = kwargs.get('file_name', 'orbit_movie.mp4')
        file_path = kwargs.get('file_path', './')
        fps = kwargs.get('fps', 25)
        keep_frames = kwargs.get('keep_frames', False)

        self.render(field_name=field_name, step=step, **kwargs)

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

        tmp_dir = kwargs.get('tmp_dir', 'temp_dir')

        if os.path.exists(tmp_dir):
            print("Error: Directory '" + tmp_dir + "' already exists. Exiting.")
            exit()

        os.mkdir(tmp_dir)

        file_name = kwargs.get('file_name', 'beltrami_instability.mp4')
        file_path = kwargs.get('file_path', './')
        fps = kwargs.get('fps', 25)
        keep_frames = kwargs.get('keep_frames', False)


        for i in range(beg, end+1):
            self.render(field_name=field_name, step=i, **kwargs)
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
        # destroy programmable filter instances
        Delete(self._prog_filter1)
        del self._prog_filter1
        Delete(self._prog_filter2)
        del self._prog_filter2

        # destroy contour instance
        Delete(self._contour)
        del self._contour

        # destroy NetCDF reader instance
        Delete(self._pvnc)
        del self._pvnc

        # update animation scene based on data timesteps
        self._animation_scene.UpdateAnimationUsingDataTimeSteps()

        ResetSession()

    def close(self):
        self._ncreader.close()
        self._clear()

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

        if field_name in self._derived_fields:
            self._contour = Contour(registrationName='Contour1', Input=self._prog_filters[field_name])
        else:
            self._contour = Contour(registrationName='Contour1', Input=self._pvnc)

        self._contour.ContourBy = ['POINTS', field_name]

        self._contour.Isosurfaces = np.linspace(vmin, vmax, n_iso, endpoint=True)
            
        self._contour.PointMergeMethod = 'Uniform Binning'

        # set active source
        SetActiveSource(self._contour)

        # show data in view
        self._contour_display = Show(self._contour, self._render_view, 'GeometryRepresentation')

        # set scalar coloring
        ColorBy(self._contour_display, ('POINTS', field_name))

        # rescale color and/or opacity maps used to include current data range
        self._contour_display.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        self._contour_display.SetScalarBarVisibility(self._render_view, True)

        # get color transfer function/color map for field_name
        self._lut = GetColorTransferFunction(field_name)

        if self._enable_opacity:
            self._lut.EnableOpacityMapping = 1
        else:
            self._lut.EnableOpacityMapping = 0

        self._lut.ApplyPreset(self.colormap, True)

        # get opacity transfer function/opacity map for field_name
        self._pwf = GetOpacityTransferFunction(field_name)

        ## Rescale transfer function
        self._lut.RescaleTransferFunction(vmin, vmax)

        if self._invert_colormap:
            self._lut.InvertTransferFunction()

        ## Rescale transfer function
        self._pwf.RescaleTransferFunction(vmin, vmax)

        points = [vmin, self._opacity_vmin, 0.5, 0.0]
        if self._opacity_points:
            for i, p in enumerate(self._opacity_points):
                v = self._opacity_values[i]
                points = points + [vmax * p, v, 0.5, 0.0]
        points = points + [vmax, self._opacity_vmax, 0.5, 0.0]

        self._pwf.Points = points


        if not self._color_vmin is None:
            points = []
            points = [vmin] + self._color_vmin # rgb
            if self._color_points:
                for i, p in enumerate(self._color_points):
                    v = self._color_values[3*i:3*i+3] # rgb
                    points = points + [p] + v
            points = points + [vmax] + self._color_vmax # rgb
            self._lut.RGBPoints = points

        if self._use_log_scale:
            # convert to log space
            self._lut.MapControlPointsToLogSpace()
            self._lut.UseLogScale = 1

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

    def _create_color_bar(self, field_name, vmin, vmax):
        self._color_bar = GetScalarBar(self._lut, self._render_view)

        if self._add_clabel:
            self._color_bar.TitleFontFamily = 'Courier'
            self._color_bar.TitleFontSize = 30
        self._color_bar.LabelFontFamily = 'Courier'
        self._color_bar.LabelFontSize = 30

        self._color_bar.ScalarBarThickness = 20
        self._color_bar.ScalarBarLength = 0.5
        self._color_bar.WindowLocation = 'Any Location'
        self._color_bar.Position = [0.89, 0.25]

        if self._add_clabel:
            if field_name in self._field_label.keys():
                self._color_bar.Title = self._field_label[field_name]
            else:
                self._color_bar.Title = field_name

            self._color_bar.TitleJustification = 'Centered'
        else:
            self._color_bar.Title = ''

        self._color_bar.AddRangeLabels = 0
        self._color_bar.RangeLabelFormat = ''

        self._color_bar.AutomaticLabelFormat = 0
        self._color_bar.LabelFormat = '%-#6.2g'

        self._color_bar.UseCustomLabels = 1
        if self._use_log_scale:
            if vmin < 1:
                cst = -len(str(int(1.0 / vmin)))
            else:
                cst = len(str(int(vmin))) - 1
            cen = len(str(int(vmax)))
            self._color_bar.CustomLabels = np.logspace(cst, cen, self._n_color_bar_ticks,
                                                       endpoint=True)
        else:
            self._color_bar.CustomLabels = np.linspace(vmin, vmax, self._n_color_bar_ticks,
                                                       endpoint=True)
        self._render_view.Update()

    def _set_camera_position(self):
        self._render_view.CameraPosition = [-8, 4, 4]
        self._render_view.CameraViewUp = [0.0, 0.0, 1.0]
        self._render_view.CameraFocalPoint = [-0.75, 0, 0]
        #self._render_view.CameraParallelScale = 2.0
        self._render_view.Update()
