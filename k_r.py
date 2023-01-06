import matplotlib.pyplot as plt
import yt

yt.enable_plugins()

f=open("QSO_agn_005.dat","r")

ds_k = []
dt_k = []

ds_r = []
dt_r = []

times_k = []
times_r = []

jet = 'JET'
qso = 'QUASAR'

for line in f.readlines()[2:]:    
    curLine=line.strip().split(" ") 
    ty = str(curLine[9])
    
    if (ty == jet):
        ds_k.append(float(curLine[5]))
        dt_k.append(float(curLine[1]))
        t=float(curLine[0])
        times_k.append(t/31557600000000000)

    if (ty == qso):
        ds_r.append(float(curLine[5]))
        dt_r.append(float(curLine[1]))
        t=float(curLine[0])
        times_r.append(t/31557600000000000)

P_k=[a/b for a,b in zip(ds_k, dt_k)]
P_r=[a/b for a,b in zip(ds_r, dt_r)]

plt.semilogy(times_k, P_k, '-', markersize=1.0, label='Kinetic')
plt.semilogy(times_r, P_r, '-', markersize=1.0, label='Radiative')
plt.ylim(1e41, 1e47)
plt.xlim(0.0,1.0)
plt.ylabel("Luminosities (erg/s)")
plt.xlabel("Time (Gyr)")
plt.savefig("k_r.png" ,dpi=300)
plt.show()
f.close()
