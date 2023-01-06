import h5py
import yt
from matplotlib import rc_context
from matplotlib.animation import FuncAnimation

ts = yt.load("/data/lizli/qso/0_1_0_05/QSO_hdf5_plt_cnt_0*")
plot = yt.SlicePlot(ts[0], "x", ("gas", "density"))
plot.annotate_timestamp(corner='upper_left', time_format='t = {time:.2f} {units}', time_unit="Myr")
plot.set_zlim(("gas", "density"), 1e-27, 2e-22)
plot.set_cmap(("gas", "density"), "algae")

fig = plot.plots[("gas", "density")].figure

def animate(i):
    ds = ts[i]
    plot._switch_ds(ds)

animation = FuncAnimation(fig, animate, frames=len(ts))

with rc_context({"mathtext.fontset": "stix"}):
    animation.save("density_x_animation_qso_0_1_0_05.mp4")
