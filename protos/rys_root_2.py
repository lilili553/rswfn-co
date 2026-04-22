# 此脚本用于找到根和权重

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

alps=np.loadtxt("alps.txt")
bets=np.loadtxt("betl.txt")
print(bets.shape)
N=10

roots=[]
weits=[]

for a,alp in enumerate(alps):
    wfunc = lambda x: np.exp(-alp*x**2)
    mu=quad(wfunc,-1,1)[0]
    Jmat=np.zeros((N,N))
    for i in range(1,N):
        Jmat[i,i-1]=np.sqrt(bets[i,a])
        Jmat[i-1,i]=np.sqrt(bets[i,a])
    
    vals,vecs=np.linalg.eigh(Jmat)
    roots.append(vals)
    weits.append(vecs[0,:]**2*mu)

roots=np.array(roots)[:,N//2:]
weits=np.array(weits)[:,N//2:]

print(roots.shape)
print(weits.shape)

np.savetxt(f"RYS_ROOTS_{N}.txt",roots.T,delimiter=", ")
np.savetxt(f"RYS_WEITS_{N}.txt",weits.T,delimiter=", ")

fig,axs=plt.subplots(1,2,figsize=(12,4))

for i in range(N//2):
    axs[0].plot(alps,roots[:,i])
    axs[0].scatter(alps,roots[:,i],s=5)
    axs[1].plot(alps,weits[:,i])
    axs[1].scatter(alps,weits[:,i],s=5)
axs[0].set_ylim(0,1)
plt.savefig(f"RYS_{N}.png",dpi=300)