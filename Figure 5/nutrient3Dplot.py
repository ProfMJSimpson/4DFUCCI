import plotly.graph_objects as go
import numpy as np
import scipy.io
import matplotlib
from matplotlib import cm
from plotly.subplots import make_subplots
import plotly.express as px

L = 2000
I = 101

cbar_custom_mat = scipy.io.loadmat('cMap',mdict=None,appendmat=True)
cbar_custom = np.array(cbar_custom_mat['cMap'], dtype=np.float32, order='K')

c_rgb = np.zeros([256,3])
for i in range(0,255):
    cbar_row = cbar_custom[i,:]
    crow = [cbar_row[0],cbar_row[1],cbar_row[2]]
    c_rgb[i,:] = crow

def mpl_to_plotly(cmap,plotly_entries):
    h = 1.0/(plotly_entries-1)
    plotly_colorscale = []

    for q in range(plotly_entries):
        C = list(map(np.uint8, cmap[q,:]*255))
        plotly_colorscale.append([q*h, 'rgb'+str((C[0],C[1],C[2]))])

    return plotly_colorscale

cool = mpl_to_plotly(c_rgb,255)

def get_slice(x,y,z, surfacecolor):
    return go.Surface(x=x,
                      y=y,
                      z=z,
                      surfacecolor=surfacecolor,
                      opacity=0.8,
                      cmin=0,
                      cmax=1,
                      colorscale=cool
                      )

# Import the workspace, use change c_x to the profile wanted. Recommended to do
# profiles one at a time to avoid issues with Plotly api
profile = scipy.io.loadmat('Prof3D\c_0',mdict=None,appendmat=True) 

# Transfer the profiles from the .mat files to python variables
c_prof = np.array(profile['c_temp'], dtype=np.float32, order='K')

xv = np.linspace(-L,L,I)
yv = np.linspace(-L,L,I)
zv = np.linspace(-L,L,I)
Xm, Ym, Zm = np.meshgrid(xv,yv,zv,indexing='ij')

##=====================================================================##
mid = int((I-1)/2 + 1)

Xm_1 = Xm[0:mid,0:mid,0:mid]
Ym_1 = Ym[0:mid,0:mid,0:mid]
Zm_1 = Zm[0:mid,0:mid,0:mid]
c_prof_1 = c_prof[0:mid,0:mid,0:mid]

Xm_2 = Xm[0:mid,0:mid,mid-1:I]
Ym_2 = Ym[0:mid,0:mid,mid-1:I]
Zm_2 = Zm[0:mid,0:mid,mid-1:I]
c_prof_2 = c_prof[0:mid,0:mid,mid-1:I]

Xm_3 = Xm[0:mid,mid-1:I,mid-1:I]
Ym_3 = Ym[0:mid,mid-1:I,mid-1:I]
Zm_3 = Zm[0:mid,mid-1:I,mid-1:I]
c_prof_3 = c_prof[0:mid,mid-1:I,mid-1:I]

Xm_4 = Xm[0:mid,mid-1:I,0:mid]
Ym_4 = Ym[0:mid,mid-1:I,0:mid]
Zm_4 = Zm[0:mid,mid-1:I,0:mid]
c_prof_4 = c_prof[0:mid,mid-1:I,0:mid]

Xm_5 = Xm[mid-1:I,0:mid,0:mid]
Ym_5 = Ym[mid-1:I,0:mid,0:mid]
Zm_5 = Zm[mid-1:I,0:mid,0:mid]
c_prof_5 = c_prof[mid-1:I,0:mid,0:mid]

Xm_6 = Xm[mid-1:I,0:mid,mid-1:I]
Ym_6 = Ym[mid-1:I,0:mid,mid-1:I]
Zm_6 = Zm[mid-1:I,0:mid,mid-1:I]
c_prof_6 = c_prof[mid-1:I,0:mid,mid-1:I]

##Xm_7 = Xm[mid-1:I,mid-1:I,mid-1:I]
##Ym_7 = Ym[mid-1:I,mid-1:I,mid-1:I]
##Zm_7 = Zm[mid-1:I,mid-1:I,mid-1:I]
##c_prof_7 = c_prof[mid-1:I,mid-1:I,mid-1:I]

Xm_8 = Xm[mid-1:I,mid-1:I,0:mid]
Ym_8 = Ym[mid-1:I,mid-1:I,0:mid]
Zm_8 = Zm[mid-1:I,mid-1:I,0:mid]
c_prof_8 = c_prof[mid-1:I,mid-1:I,0:mid]

opac = [[0, 1],[0.01, 1],[0.2,1],[0.25,0.4],[0.4,0.3],[0.8, 0.6],[1, 0.6]]

fig5a = go.Volume(
    x=Xm_1.flatten(),
    y=Ym_1.flatten(),
    z=Zm_1.flatten(),
    value=c_prof_1.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5b = go.Volume(
    x=Xm_2.flatten(),
    y=Ym_2.flatten(),
    z=Zm_2.flatten(),
    value=c_prof_2.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5c = go.Volume(
    x=Xm_3.flatten(),
    y=Ym_3.flatten(),
    z=Zm_3.flatten(),
    value=c_prof_3.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5d = go.Volume(
    x=Xm_4.flatten(),
    y=Ym_4.flatten(),
    z=Zm_4.flatten(),
    value=c_prof_4.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5e = go.Volume(
    x=Xm_5.flatten(),
    y=Ym_5.flatten(),
    z=Zm_5.flatten(),
    value=c_prof_5.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5f = go.Volume(
    x=Xm_6.flatten(),
    y=Ym_6.flatten(),
    z=Zm_6.flatten(),
    value=c_prof_6.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

fig5g = go.Volume(
    x=Xm_8.flatten(),
    y=Ym_8.flatten(),
    z=Zm_8.flatten(),
    value=c_prof_8.flatten(),
    isomin=0,
    isomax=1,
    opacity=0.3, # needs to be small to see through all surfaces
    surface_count=70, # needs to be a large number for good volume rendering
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorscale=cool
    )

x2d = np.linspace(0,L,mid)
y2d = np.linspace(0,L,mid)
z2d = np.linspace(0,L,mid)

## XY plane
Xmxy,Ymxy = np.meshgrid(x2d,y2d)
Zmxy = np.zeros(Xmxy.shape)
surfcolorxy = c_prof[mid-1:I,mid-1:I,mid-1]
xyslice = get_slice(Xmxy,Ymxy,Zmxy,surfcolorxy)

## YZ plane
Ymyz,Zmyz = np.meshgrid(y2d,z2d)
Xmyz = np.zeros(Ymyz.shape)
surfcoloryz = c_prof[mid-1,mid-1:I,mid-1:I]
yzslice = get_slice(Xmyz,Ymyz,Zmyz,surfcoloryz)

## XZ plane
Xmxz,Zmxz = np.meshgrid(x2d,z2d)
Ymxz = np.zeros(Xmxz.shape)
surfcolorxz = c_prof[mid-1:I,mid-1,mid-1:I]
xzslice = get_slice(Xmxz,Ymxz,Zmxz,surfcolorxz)


fig5h = xyslice
fig5i = yzslice
fig5j = xzslice

fig5 = make_subplots(specs=[[{"secondary_y": False}]])
fig5.add_trace(fig5a)
fig5.add_trace(fig5b)
fig5.add_trace(fig5c)
fig5.add_trace(fig5d)
fig5.add_trace(fig5e)
fig5.add_trace(fig5f)
fig5.add_trace(fig5g)
fig5.add_trace(fig5h)
fig5.add_trace(fig5i)
fig5.add_trace(fig5j)

fig5.update_layout(scene_xaxis_showticklabels=False,
                  scene_yaxis_showticklabels=False,
                  scene_zaxis_showticklabels=False,
                  scene_aspectmode='manual',
                  scene_aspectratio=dict(x=1,y=1,z=1),
                  paper_bgcolor='rgba(0,0,0,0)',
                  plot_bgcolor='rgba(0,0,0,0)',
                  scene = dict(
                    xaxis = dict(
                         backgroundcolor="rgb(255,255,255)",
                         gridcolor="Grey",
                         showbackground=True,
                         zerolinecolor="white",
                         range = [-0.2*L,0.2*L],
                         tickmode = 'linear',
                         tick0 = -400,
                         dtick = 200),
                    yaxis = dict(
                        backgroundcolor="rgb(255,255,255)",
                        gridcolor="Grey",
                        showbackground=True,
                        zerolinecolor="white",
                        range = [-0.2*L,0.2*L],
                        tickmode = 'linear',
                        tick0 = -400,
                        dtick = 200),
                    zaxis = dict(
                        backgroundcolor="rgb(255,255,255)",
                        gridcolor="Grey",
                        showbackground=True,
                        zerolinecolor="white",
                        range = [-0.2*L,0.2*L],
                        tickmode = 'linear',
                        tick0 = -400,
                        dtick = 200),
                    xaxis_title='',
                    yaxis_title='',
                    zaxis_title=''),
                  scene_camera = dict(
                        up=dict(x=0, y=0, z=1),
                        center=dict(x=0, y=0, z=0),
                        eye=dict(x=1.261249796,y=1.503098975,z=0.9149685367) #dict(x=1.1414, y=1.4875, z=0.7825)
                        ))
fig5.layout.scene.camera.projection.type = "orthographic"
fig5.show()

# When the figure appears, save as a PNG into the Nutrient3D page
