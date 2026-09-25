import math, cmath, numpy as np
from matplotlib.colors import hsv_to_rgb
AREF=0.45
def colour(z,aref=AREF):
    ph=(cmath.phase(z)/(2*math.pi))%1.0; a=min(abs(z)/aref,1.0)
    rgb=hsv_to_rgb((ph,0.85,0.95)); return tuple(1-(1-rgb)*a)
def wheel(ax,fs=6.5):
    th=np.linspace(0,2*np.pi,181); r=np.linspace(0,1,25); T,Rr=np.meshgrid(th,r)
    C=np.zeros(T.shape+(3,))
    for i in range(T.shape[0]):
        for j in range(T.shape[1]): C[i,j]=colour(Rr[i,j]*AREF*cmath.exp(1j*T[i,j]))
    ax.pcolormesh(T,Rr,C[:-1,:-1],shading='flat',rasterized=True); ax.set_yticks([])
    ax.set_xticks([k*np.pi/4 for k in range(0,8,2)]); ax.set_xticklabels(['0','$\\pi/2$','$\\pi$','$3\\pi/2$'],fontsize=fs)
    ax.tick_params(pad=-2); ax.spines['polar'].set_linewidth(.5)
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
