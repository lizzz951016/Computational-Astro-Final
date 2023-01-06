import h5py
import yt
import numpy as np
import matplotlib.pyplot as plt

cool_gas_005 = []
labels_005 = []

cool_gas_05 = []
labels_05 = []

cool_gas_1 = []
labels_1 = []

ts_005 = yt.load("/data/lizli/qso/test/QSO_hdf5_plt_cnt_0*")

for ds in ts_005:
    ad = ds.all_data()
    cool_ad = ad.cut_region(['obj["gas", "temperature"] < 5e5'])
#    cool_gas_k.append(cool_ad["cell_mass"].in_units("Msun").sum())
    cool_gas_005.append(cool_ad.quantities.total_quantity("cell_mass").in_units("Msun"))
    labels_005.append(float("%f" % ds.current_time.in_units("Gyr")))


ts_05 = yt.load("/data/lizli/qso/0_5/QSO_hdf5_plt_cnt_0*")

for ds in ts_05:
    ad = ds.all_data()
    cool_ad = ad.cut_region(['obj["gas", "temperature"] < 5e5'])
#    cool_gas_k.append(cool_ad["cell_mass"].in_units("Msun").sum())
    cool_gas_05.append(cool_ad.quantities.total_quantity("cell_mass").in_units("Msun"))
    labels_05.append(float("%f" % ds.current_time.in_units("Gyr")))


ts_1 = yt.load("/data/lizli/qso/1/QSO_hdf5_plt_cnt_0*")

for ds in ts_1:
    ad = ds.all_data()
    cool_ad = ad.cut_region(['obj["gas", "temperature"] < 5e5'])
#    cool_gas_k.append(cool_ad["cell_mass"].in_units("Msun").sum())
    cool_gas_1.append(cool_ad.quantities.total_quantity("cell_mass").in_units("Msun"))
    labels_1.append(float("%f" % ds.current_time.in_units("Gyr")))



plt.semilogy(labels_005, cool_gas_005, '-',  label="$\epsilon_{r}$=0.05")
plt.semilogy(labels_05, cool_gas_05, '-',  label="$\epsilon_{r}$=0.5")
plt.semilogy(labels_1, cool_gas_1, '-',  label="$\epsilon_{r}$=1")

plt.ylim(1e10, 1e12)
plt.xlim(0.0,1.0)
plt.ylabel("M$_{cold}$ (M$_{\odot}$)", fontsize=15)
plt.xlabel("Time (Gyr)", fontsize=15)
plt.title("Cold Gas")
plt.legend(loc="upper left")
plt.savefig("cold_gas.png" ,dpi=300)
plt.show()
