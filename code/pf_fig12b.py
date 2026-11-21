"""Dipolar row of Fig. 12: replica means with standard errors (dip_binder.py), Binder cumulant, single-state phase fields."""
from pubstyle import *
import pickle, glob, math, cmath, numpy as np, os, io, contextlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.ticker as mt
from pcolor import colour, outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
L=48
if os.path.exists('dip/snaps48b.pkl'): snaps=pickle.load(open('dip/snaps48b.pkl','rb'))
else:
    from dip_mc import simulate
    with contextlib.redirect_stdout(io.StringIO()):
        snaps=simulate(L,'point',[1e9,0.6,0.45,0.3],3000,500,seed=21,tag='snap')
    snaps=[dict(T=s['T'],h=s['h']) for s in snaps]; pickle.dump(snaps,open('dip/snaps48b.pkl','wb'))
R={}
for f in glob.glob('dip/binder_*_s*.pkl'):
    Lx=int(f.split('_')[1]); R.setdefault(Lx,[]).append(pickle.load(open(f,'rb')))
Ls=sorted(R); Ts=np.array([r['T'] for r in R[Ls[0]][0]]); Tp=np.where(Ts>1e3,4.0,Ts)
cols={24:'#c6dbef',32:'#6baed6',48:'#2171b5',64:'#08306b'}
def ms(Lx,key):
    v=np.array([[run[i][key] for i in range(len(Ts))] for run in R[Lx]]); return v.mean(0), v.std(0,ddof=1)/math.sqrt(len(v))
def tax(a):
    a.set_xscale('log'); a.set_xticks([0.4,0.6,1,4]); a.set_xticklabels(['0.4','0.6','1','$\\infty$']); a.xaxis.set_minor_formatter(mt.NullFormatter()); a.set_xlabel('$T/D$',labelpad=1)
    a.axvspan(0.43,0.47,color='#dddddd',lw=0,zorder=0)
fig=plt.figure(figsize=(TWO,1.75))
gs=fig.add_gridspec(1,4,wspace=0.55,left=0.07,right=0.995,top=0.86,bottom=0.22)
a=fig.add_subplot(gs[0,0])
for Lx in Ls:
    m,e=ms(Lx,'vabs'); a.errorbar(Tp,m,yerr=e,fmt='o-',color=cols[Lx],ms=2,lw=.8,elinewidth=.6,capsize=1,label=f'$L={Lx}$')
tax(a); a.set_ylabel('$\\overline{|\\langle v\\rangle|}$ (center)'); a.legend(fontsize=6,handlelength=1,loc='upper right'); a.set_title('(d)',loc='left',pad=2)
a=fig.add_subplot(gs[0,1])
for Lx in Ls:
    m,e=ms(Lx,'U'); a.errorbar(Tp,m,yerr=e,fmt='o-',color=cols[Lx],ms=2,lw=.8,elinewidth=.6,capsize=1)
a.axhline(0,color='k',lw=.4); a.axhline(0.5,color='#888',ls=':',lw=.8)
tax(a); a.set_ylabel('Binder cumulant $U$',labelpad=1); a.set_title('(e)',loc='left',pad=2)
for j,k in enumerate([0,3]):
    a=fig.add_subplot(gs[0,2+j]); h=snaps[k]['h']; P0={}
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
    a.set_title('('+'fg'[j]+') $\\arg\\Psi$, '+('$T=\\infty$' if T>1e3 else f'$T={T}D$'),loc='left',pad=2)
save(fig,'fig12_dipolar_row')
