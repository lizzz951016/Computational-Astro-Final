import matplotlib.pyplot as plt
import yt

yt.enable_plugins()

f=open("QSO_agn_005.dat","r")

ds = []
dt = []
cooling = []
labels = []
times = []
h = [0]
s = [0]
p = [0]
average = []
cooling_rate = []

energy=0
i=(3155760000000000/2)
j=(3155760000000000/2)
k=0

for line in f.readlines()[2:]:    
    curLine=line.strip().split(" ")
    energy+=float(curLine[5])   
    ds.append(float(curLine[5]))
    dt.append(float(curLine[1]))
    t=float(curLine[0])
    times.append(t/31557600000000000)
    if t>j:
        p.append(t)
        s.append(((t/(31557600000000000))-h[k])/2+h[k])
        h.append(t/(31557600000000000))
        energy=energy/(t-p[k])
        if s[0]==0:
            del s[0]
            del h[0]
            del p[0]
            k-=1
        average.append(energy)
        k+=1
        j+=i
        energy=0

P=[a/b for a,b in zip(ds, dt)]

ts = yt.load("/data/lizli/qso/test/QSO_hdf5_plt_cnt_0*")

for ds in ts:
    sp = ds.sphere("c", (100, "kpc"))
    cooling.append(sp.quantities.total_quantity("xray_luminosity"))
    labels.append(float("%f" % ds.current_time.in_units("Gyr")))

plt.semilogy(labels, cooling, '--', color='black', label='Cooling Rate')
plt.semilogy(times, P, 'o', color='black', markersize=0.5, label='AGN power')
plt.semilogy(s[5:], average[5:], '-s', color='red', markersize=3, linewidth=1,  label='Averaged total power')
plt.ylim(5e43, 5e47)
plt.xlim(0.0,1.0)
plt.ylabel("Luminosities (erg/s)")
plt.xlabel("Time (Gyr)")
plt.legend()
plt.savefig("Figure_1.png" ,dpi=300)
plt.show()
f.close()
