import pickle, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from square_mc import setup, vfield
out=pickle.load(open(sys.argv[1],'rb')); L=int(sys.argv[2]); idx=[int(x) for x in sys.argv[3].split(',')]
X,Y,sgn,R=setup(L,0,0)
fig,ax=plt.subplots(1,len(idx),figsize=(6.5*len(idx),6.5))
for a,i in zip(ax,idx):
    h=out[i]['h']; v,Hd,Vd=vfield(h,X,Y,sgn,L)
    for x in range(L):
        for y in range(L):
            if Hd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),1.84,.84,fc='#e8743a' if v[x,y]>0 else '#2f6db5',ec='k',lw=.4))
            if Vd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),.84,1.84,fc='#e8743a' if v[x,y]>0 else '#2f6db5',ec='k',lw=.4))
    a.set_xlim(L/2-10,L/2+10); a.set_ylim(L/2-10,L/2+10); a.set_aspect('equal'); a.set_title(f"T = {out[i]['T']:g} (chessboard, centre)")
plt.tight_layout(); plt.savefig(sys.argv[4],dpi=100)
# structure of domino occupancy at the last index
h=out[idx[-1]]['h']; v,Hd,Vd=vfield(h,X,Y,sgn,L)
for name,A in (('H',Hd),('V',Vd)):
    F=np.abs(np.fft.fft2(A.astype(float)-A.mean()))**2; F[0,0]=0
    top=np.argsort(F.ravel())[::-1][:6]; print(name,[(round(t//L*2/L,3),round(t%L*2/L,3),round(F.ravel()[t]/F.sum(),3)) for t in top],'(q/π)')
