import pickle, sys, numpy as np, math, cmath
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from square_mc import setup, vfield
from core import rect, short_cell
from plot_defect_maps import outline
from plot_psi import colour
out=pickle.load(open(sys.argv[1],'rb')); L=int(sys.argv[2]); idx=[int(x) for x in sys.argv[3].split(',')]
X,Y,sgn,R=setup(L,0,0)
fig,ax=plt.subplots(1,len(idx),figsize=(7*len(idx),7))
for a,i in zip(np.atleast_1d(ax),idx):
    o=out[i]; h=o['h']; v,Hd,Vd=vfield(h,X,Y,sgn,L)
    for x in range(L):
        for y in range(L):
            ds=[]
            if Hd[x,y]: ds.append(((x,y),(x+1,y)))
            if Vd[x,y]: ds.append(((x,y),(x,y+1)))
            for d in ds:
                s=short_cell(*d); l=d[1] if s==d[0] else d[0]; vv=v[x,y]
                a.add_patch(Polygon(outline(list(rect(s))),closed=True,fc='#e8743a' if vv>0 else '#2f6db5',ec='none'))
                a.add_patch(Polygon(outline(list(rect(l))),closed=True,fc='#f6c6a8' if vv>0 else '#aac6ea',ec='none'))
                a.add_patch(Polygon(outline(list(rect(s))+list(rect(l))),closed=True,fill=False,ec='k',lw=.4))
    a.set_aspect('equal'); a.autoscale_view(); a.axis('off'); a.set_title(f"T = {o['T']:g}: E/L² = {o['E']:.4f}",fontsize=12)
plt.tight_layout(); plt.savefig(sys.argv[4],dpi=110); print('ok')
