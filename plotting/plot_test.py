## mayavi2 -x volume_rendering.py -o
import numpy as np
import matplotlib.pyplot as plt
from tools.nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid

mpl.rcParams['font.size'] = 16

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.3f"):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')

    label = r"t = " + fmt % (time)

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox
    )

grid = '32'
step = 62
opacity = 0.01
surface_count = 150

ncreader = nc_reader()
ncreader.open('../examples/beltrami_' + grid + '_fields.nc')

t = ncreader.get_all('t')

origin = ncreader.get_box_origin()
extent = ncreader.get_box_extent()
ncells = ncreader.get_box_ncells()

xmin = origin[0]
xmax = origin[0] + extent[0]
ymin = origin[1]
ymax = origin[1] + extent[1]
zmin = origin[2]
zmax = origin[2] + extent[2]

##pres = ncreader.get_dataset(step=step, name='pressure')
x_vor = ncreader.get_dataset(step=step, name='x_vorticity')
y_vor =	ncreader.get_dataset(step=step, name='y_vorticity')
z_vor =	ncreader.get_dataset(step=step, name='z_vorticity')

#nz, ny, nx = x_vor.shape

vor = np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)

#conda install -c conda-forge plotly
#conda install -c conda-forge python-kaleido



import plotly.graph_objects as go
import numpy as np

# 12 July 2022
# https://plotly.com/python/renderers/
#import plotly.io as pio
#pio.renderers.default = "png"

# 12 July 2022
# https://plotly.com/python/v3/matplotlib-colorscales/#formatting-the-colormap
def matplotlib_to_plotly(cmap, pl_entries):
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        #C = map(np.uint8, np.array(cmap(k*h)[:3])*255)
        C = np.array(cmap(k*h)[:3])*255
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

    return pl_colorscale

kbc_cmap = cc.cm['kbc']

kbc = matplotlib_to_plotly(kbc_cmap, 256)


#help(go.Volume)
#exit()

def copy_periodic_layers(field):
        nz, ny, nx = field.shape
        field_copy = np.empty((nz, ny+1, nx+1))
        field_copy[:, 0:ny, 0:nx] = field.copy()
        field_copy[:, ny, :] = field_copy[:, 0, :]
        field_copy[:, :, nx] = field_copy[:, :, 0]
        return field_copy


x = ncreader.get_all(name='x')
y = ncreader.get_all(name='y')
z = ncreader.get_all(name='z')
vor = copy_periodic_layers(vor)

# add periodic values
x = np.append(x, abs(x[0]))
y = np.append(y, abs(y[0]))

#vor[:, :, :] = 0.0
#vor[:, 0, 0] = y

vor = np.transpose(vor, axes=[2, 1, 0])

#vor[:, :, :] = 0.0
#vor[0, 0, :] = y

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# 13 July 2022
# https://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d
assert np.all(X[:,0,0] == x)
assert np.all(Y[0,:,0] == y)
assert np.all(Z[0,0,:] == z)



fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=vor.flatten(),
    #isomin=0.01,
    #isomax=0.99,
    opacity=opacity, # needs to be small to see through all surfaces
    surface_count=surface_count, # needs to be a large number for good volume rendering
    colorscale=kbc,
    colorbar=dict(
            #title=dict(text='$|\omega|$', side='top'),
            thickness = 20,
            len=0.75,
            xpad=5,
            orientation='v',
            tickfont=dict(
                size=22,
                color="black"
            ),
        )
    ))


tickvals = np.array([-1.5, -0.75, 0.0, 0.75, 1.5])
ticktext = [' -3/2 ', ' -3/4 ', ' 0 ', ' 3/4 ', ' 3/2 ']

# 12 July 2022
# https://plotly.com/python/3d-camera-controls/
camera = dict(
    up=dict(x=0, y=0, z=1),
    center=dict(x=0, y=0, z=-0.1),
    eye=dict(x=1.4, y=1.4, z=1.4)
)


fig.update_layout(
    scene_camera=camera,
    scene = dict(
                    xaxis = dict(
                        tickmode = 'array',
                        tickvals = tickvals,
                        ticktext = ticktext
                    ),
                    yaxis = dict(
                        tickmode = 'array',
                        tickvals = tickvals,
                        ticktext = ticktext
                    ),
                    zaxis = dict(
                        tickmode = 'array',
                        tickvals = tickvals,
                        ticktext = ticktext
                    ),
                    xaxis_title=dict(text=r'x', font=dict(size=24, color='black')),
                    yaxis_title=dict(text=r'y', font=dict(size=24, color='black')),
                    zaxis_title=dict(text=r'z', font=dict(size=24, color='black'))),
                    margin=dict(r=10, b=5, l=10, t=5),
    font=dict(
        family="Arial", # Times New Roman
        size=16,
        color="black"
    ),
)


fig.write_image("test_plotly.png", scale=1.5, width=1024, height=1024)

image = plt.imread("test_plotly.png", format='png')
plt.figure(figsize=(9, 9))
plt.imshow(image)
plt.axis('off')
add_timestamp(plt, t[step], xy=(0.05, 0.95), fmt="%.2f")

plt.tight_layout()
plt.savefig('test_plotly.eps', format='eps')

ncreader.close()
