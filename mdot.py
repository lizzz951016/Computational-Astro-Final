import matplotlib.pyplot as plt
import yt

yt.enable_plugins()

f=open("QSO_agn_005.dat","r")

mdot_005 = []
times_005 = []

for line in f.readlines()[2:]:    
    curLine=line.strip().split(" ") 
    mdot_005.append(float(curLine[3]))
    t=float(curLine[0])
    times_005.append(t/31557600000000000)

f=open("QSO_agn_05.dat","r")

mdot_05 = []
times_05 = []

for line in f.readlines()[2:]:    
    curLine=line.strip().split(" ") 
    mdot_05.append(float(curLine[3]))
    t=float(curLine[0])
    times_05.append(t/31557600000000000)

f=open("QSO_agn_1.dat","r")

mdot_1 = []
times_1 = []

for line in f.readlines()[2:]:    
    curLine=line.strip().split(" ") 
    mdot_1.append(float(curLine[3]))
    t=float(curLine[0])
    times_1.append(t/31557600000000000)

plt.semilogy(times_005, mdot_005, '-', markersize=2, label='$\epsilon_{r}$=0.05')
#plt.semilogy(times_05, mdot_05, '-', markersize=0.5, label='$\epsilon_{r}$=0.5')
#plt.semilogy(times_1, mdot_1, '-', markersize=0.5, label='$\epsilon_{r}$=1')
plt.ylim(1e-4, 1e4)
plt.xlim(0.0,1.0)
plt.ylabel("$\dot{M}_{BH}$ ($M_{\odot} yr^{-1}$)", fontsize=15)
plt.xlabel("Time (Gyr)", fontsize=15)
plt.legend()
plt.savefig("mdot_005.png" ,dpi=300)
plt.show()
f.close()
