from pubstyle import *
import pickle, math, cmath, numpy as np, os, io, contextlib
from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection
from pcolor import colour, outline
from square_mc import setup, vfield
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
L=48
if os.path.exists('dip/snaps48.pkl'): snaps=pickle.load(open('dip/snaps48.pkl','rb'))
else:
    from dip_mc import simulate
    with contextlib.redirect_stdout(io.StringIO()):
        snaps=simulate(L,'point',[1e9,0.8,0.6,0.5,0.45,0.4,0.3],3000,500,seed=21,tag='snap')
    pickle.dump([dict(T=s['T'],h=s['h']) for s in snaps],open('dip/snaps48.pkl','wb'))
R={Lx:pickle.load(open(f'dip/fss_{Lx}_point.pkl','rb')) for Lx in (24,32,48,64)}
Ts=np.array([r['T'] for r in R[24]]); Tp=np.where(Ts>1e3,4.0,Ts)
fig=plt.figure(figsize=(TWO,3.45))
gs=fig.add_gridspec(2,4,height_ratios=[1,1],wspace=0.42,hspace=0.25,left=0.07,right=0.99,top=0.95,bottom=0.02)
cols={24:'#c6dbef',32:'#6baed6',48:'#2171b5',64:'#08306b'}
import matplotlib.ticker as mt
def tax(a):
    a.set_xscale('log'); a.set_xticks([0.2,0.5,1,2,4]); a.set_xticklabels(['0.2','0.5','1','2','$\\infty$']); a.xaxis.set_minor_formatter(mt.NullFormatter()); a.set_xlabel('$T/D$',labelpad=1)
a=fig.add_subplot(gs[0,0])
for Lx in R: a.plot(Tp,[r['C'] for r in R[Lx]],'o-',color=cols[Lx],ms=2,lw=.8,label=f'$L={Lx}$')
tax(a); a.set_ylabel('specific heat'); a.legend(fontsize=6,handlelength=1,loc='upper right'); a.set_title('(a)',loc='left',pad=2)
a=fig.add_subplot(gs[0,1])
for Lx in R: a.plot(Tp,[np.mean(np.abs(r['vmap'])) for r in R[Lx]],'o-',color=cols[Lx],ms=2,lw=.8)
tax(a); a.set_ylabel('$\\overline{|\\langle v\\rangle|}$ (center)'); a.set_title('(b)',loc='left',pad=2)
a=fig.add_subplot(gs[0,2]); Ls=np.array(list(R))
x2=[-np.polyfit(np.log(Ls),np.log([R[Lx][i]['m2'] for Lx in R]),1)[0] for i in range(len(Ts))]
a.plot(Tp,x2,'o-',color=CAT[0],ms=2.5,lw=.8); a.axhline(x2[0],color='#888',ls=':',lw=.8); a.axhline(0,color='k',lw=.4)
tax(a); a.set_ylabel('effective $2x$'); a.set_title('(c)',loc='left',pad=2)
a=fig.add_subplot(gs[0,3]); X,Y,sgn,Rr=setup(L,0,0); h=snaps[-1]['h']; v,Hd,Vd=vfield(h,X,Y,sgn,L)
for x in range(L):
    for y in range(L):
        if Hd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),1.84,.84,fc=OR if v[x,y]>0 else BL,ec='k',lw=.25))
        if Vd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),.84,1.84,fc=OR if v[x,y]>0 else BL,ec='k',lw=.25))
a.set_xlim(L/2-6,L/2+6); a.set_ylim(L/2-6,L/2+6); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title('(d) crystal, $T=0.3D$',loc='left',pad=2)
for j,k in enumerate([0,2,4,6]):
    a=fig.add_subplot(gs[1,j]); h=snaps[k]['h']; P0={}
    for x in range(L):
        for y in range(L): P0[(x,y)]=0.5*sum(cmath.exp(1j*math.pi*h[x+aa,y+bb]/4) for aa in (0,1) for bb in (0,1))
    polys=[];fcs=[]
    for (x,y),z in P0.items():
        zz=0;ws=0
        for dx in range(-3,4):
            for dy in range(-3,4):
                q=(x+dx,y+dy)
                if q in P0: w=math.exp(-(dx*dx+dy*dy)/(2*1.3**2)); zz+=w*P0[q]; ws+=w
        polys.append(Polygon(outline(list(rect((x,y)))),closed=True)); fcs.append(colour(zz/ws*0.45/0.6))
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors=fcs,linewidths=.2))
    a.set_aspect('equal'); a.autoscale_view(); a.margins(0.01); a.axis('off'); T=snaps[k]['T']
    Ts_=r'\infty' if T>1e3 else str(T); a.set_title('('+'efgh'[j]+') $T='+Ts_+('$' if T>1e3 else 'D$'),loc='left',pad=1)
save(fig,'fig12_dipolar')
