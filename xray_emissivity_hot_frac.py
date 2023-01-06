import yt
import h5py
from yt.units import kpc
import numpy as np
import math
import matplotlib.pyplot as plt
from yt.visualization.base_plot_types import get_multi_plot
from matplotlib.colors import LogNorm
import matplotlib as mpl


ds = yt.load("/data/lizli/qso/0_05/QSO_hdf5_plt_cnt_0133")

orient = "vertical"

SB0_arr = [0.0 for i in range(905)]
N_arr = [0.0 for i in range(905)]
dr=0.8
f = [[0.0]*1024 for i in range(1024)]
f=np.array(f)

fig, axes, colorbars = get_multi_plot(1, 1, colorbar=orient, bw=4)


ds.periodicity=(True,True,True)
ad = ds.all_data()
hot = ad.cut_region(['(obj["temperature"] > 5.8e6)'])

proj = yt.ProjectionPlot(ds, "x", ("gas", "xray_emissivity"), data_source=hot)
proj_frb = proj.data_source.to_frb((160, "kpc"), 1024)
var_frb = np.array(proj_frb[("gas", "xray_emissivity")])

xray_axes = [axes[0][0]]

for xax in xray_axes:
    xax.xaxis.set_visible(False)
    xax.yaxis.set_visible(False)

for iy in range(1024):
    for ix in range(1024):
        r=((ix-511.5)**2+(iy-511.5)**2)**(0.5)
        ir=math.floor(r/dr)
        SB0_arr[ir]+=var_frb[iy,ix]
        N_arr[ir]+=1

for i in range(905):
    SB0_arr[i]=SB0_arr[i]/N_arr[i]

for i in range(1024):
    for j in range(1024):
        r=((i-511.5)**2+(j-511.5)**2)**(0.5)
        ir=math.floor(r/dr)
        f[i,j]=(var_frb[i,j]-SB0_arr[ir])/SB0_arr[ir]

plots = [
    xray_axes[0].imshow(f, origin="lower"),
]

plots[0].set_cmap("GREEN")

titles = ["xray emissivity (fractional variation)"]

plots[0].set_clim(-0.2,0.2)

for p, cax, t in zip(plots, colorbars, titles):
    cbar = fig.colorbar(p, cax=cax, orientation=orient)
    cbar.set_label(t)

fig.savefig("xray_emissivity_frac_1")
