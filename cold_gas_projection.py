import h5py
import yt
import numpy as np


yt.enable_plugins()

ds = yt.load("/data/lizli/qso/test/QSO_hdf5_plt_cnt_0130")

ds.periodicity=(True,True,True)
ad1 = ds.all_data()
cold = ad1.cut_region(['obj["temperature"] < 5e5'])
slc = yt.ProjectionPlot(ds, "x",  ("xray_emissivity"), width=(80.0, "kpc"), data_source=cold)
slc.annotate_timestamp(corner='upper_left', time_format='t = {time:.2f} {units}', time_unit="Gyr", text_args={'color':'white'})


slc.set_colorbar_label("xray_emissivity", "Projected Xray Emissivity")
slc.set_font({"size": 45})
slc.hide_axes(draw_frame=True)
slc.set_background_color(field='xray_emissivity', color='k')
slc.set_cmap(("gas", "xray_emissivity"), "gist_heat")
slc.save()