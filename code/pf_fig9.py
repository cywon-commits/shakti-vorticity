from pubstyle import *
import pickle, math, cmath, numpy as np, io, contextlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pcolor import outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
L=48
def load(H,k): return pickle.load(open(f'hole/s_{L}_{H}_k{k}.pkl','rb'))
prof=pickle.load(open('hole/prof.pkl','rb'))
base=load(0,0)
def ring(V,cell,d0,d1):
    a=0j;b=0j;n=0
    for i in range(L):
        for j in range(L):
            if not cell[i,j]: continue
            d=max(abs(i+.5-L/2),abs(j+.5-L/2))
            if d0<=d<d1: n+=1; a+=V[i,j]*cmath.exp(-1j*math.pi/2*(i+j)); b+=V[i,j]*cmath.exp(-1j*math.pi/2*(i-j))
    return a/n,b/n
fig=plt.figure(figsize=(TWO,3.6))
gs=fig.add_gridspec(2,3,width_ratios=[1.35,1,1],wspace=0.06,hspace=0.12,left=0.075,right=0.93,top=0.95,bottom=0.1)
ap=fig.add_subplot(gs[0,0]); aa=fig.add_subplot(gs[1,0],sharex=ap)
style={(8,0):(CAT[0],'o','$H=8$, $k=0$'),(8,1):(CAT[1],'s','$H=8$, $k=+1$'),(8,-1):(CAT[1],'s',None),(10,0):(CAT[3],'^','$H=10$, $k=0$'),(10,1):(CAT[2],'D','$H=10$, $k=+1$')}
for key_,rows in prof.items():
    c,m,lab=style[key_]
    ap.plot([r[0] for r in rows],[r[1] for r in rows],'-',marker=m,color=c,ms=2.5,lw=.8,label=lab)
    aa.plot([r[0] for r in rows],[r[3] for r in rows],'-',marker=m,color=c,ms=2.5,lw=.8)
R0=[(d0+1,*ring(base['vmean'],base['cell'],d0,d0+2)) for d0 in np.arange(0,L/2,2)]
ap.plot([r[0] for r in R0],[math.sqrt(2*(abs(r[1])**2+abs(r[2])**2)) for r in R0],'k--',lw=.8,label='no hole')
ap.set_ylabel('local amplitude'); plt.setp(ap.get_xticklabels(),visible=False)
ap.set_ylim(0.35,0.8)
h_,l_=ap.get_legend_handles_labels(); aa.legend(h_,l_,fontsize=6.5,loc='lower right',handlelength=1.4)
aa.axhline(0,color='#999',lw=.5); aa.set_ylim(-1.1,1.1); aa.set_xlabel('distance from center (cells)'); aa.set_ylabel('alignment $\\cos\\Delta\\varphi$')
ap.set_title('(a)',loc='left',pad=2); aa.set_title('(b)',loc='left',pad=2)
nv=Normalize(-.6,.6)
panels=[((0,0),'(c) no hole'),((8,0),'(d) $H=8$, $k=0$'),((8,1),'(e) $H=8$, $k=+1$'),((10,0),'(f) $H=10$, $k=0$')]
for idx,((key_,title)) in enumerate(panels):
    D=load(*key_); a=fig.add_subplot(gs[idx//2,1+idx%2]); V=D['vmean']; W=13; polys=[]; fcs=[]; xs=[];ys=[]
    for i in range(L):
        for j in range(L):
            if D['cell'][i,j] and abs(i+.5-L/2)<=W and abs(j+.5-L/2)<=W:
                r=list(rect((i,j))); polys.append(Polygon(outline(r),closed=True)); fcs.append(CMAPV(nv(V[i,j]))); xs+=[p[0] for p in r]; ys+=[p[1] for p in r]
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors=fcs,linewidths=.2))
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off'); a.set_title(title,loc='left',pad=2)
    last=a
cax=last.inset_axes([1.0,0.62,0.045,0.36]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=CMAPV),cax=cax,ticks=[-.6,0,.6])
cb.ax.tick_params(labelsize=6,length=2,width=.5,pad=1); cb.outline.set_linewidth(.5); cb.set_label('$\\langle v\\rangle$',fontsize=7,labelpad=1)
save(fig,'fig9_holes')
