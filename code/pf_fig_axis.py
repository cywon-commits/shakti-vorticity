"""Straight-edged Shakti samples (data: open_vortex.py, axis_psi.py 48)."""
from pubstyle import *
import pickle, math, cmath, io, contextlib, numpy as np
from matplotlib.patches import Polygon, Rectangle, FancyArrow
from matplotlib.collections import PatchCollection
from pcolor import colour, wheel, outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, pcen
from axis_region import CELLS
D=pickle.load(open('exact/axis_psi_48.pkl','rb')); G=pickle.load(open('exact/open_gs.pkl','rb')); F=pickle.load(open('exact/open_forced.pkl','rb'))
G[(5,'forced')]=dict(box=F[5]['box'],A=F[5]['box'],bbox=F[5]['bbox'],**{k:F[5]['forced'][k] for k in ('isl','spins','labs')})
SAT=1.3066; amp_norm=Normalize(0,1); EXTRA='#f2b48c'
black=lambda c:(c[0]+c[1])%2==0
def draw(a,key,mode):
    d=D[key]; P=d['P']; polys=[]; fcs=[]
    for c in d['box']:
        polys.append(Polygon(outline(list(CELLS[c])),closed=True))
        if c not in P: fcs.append(EXTRA)
        else: fcs.append(plt.cm.viridis(amp_norm(abs(P[c])/SAT)) if mode=='amp' else colour(P[c]/max(abs(P[c]),1e-9)*0.45))
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors=fcs,linewidths=.3))
    x0,y0=d['origin']; X=d['X']
    a.set_xlim(x0-1,x0+X+1); a.set_ylim(y0-1,y0+X+1); a.set_aspect('equal'); a.axis('off')
def gs_panel(a,key):
    g=G[(5,key)]; x0,y0,x1,y1=g['bbox']
    for c in g['box']:
        a.add_patch(Polygon(outline(list(CELLS[c])),closed=True,fc=(EXTRA if c not in g['A'] else ('#b9c6de' if black(c) else '#dfe6f2')),ec='none'))
    for k,(p,q) in enumerate(g['isl']):
        s=g['spins'][k]; a0,b0=(p,q) if s==1 else (q,p)
        dx,dy=b0[0]-a0[0],b0[1]-a0[1]
        a.add_patch(FancyArrow(a0[0]+.2*dx,a0[1]+.2*dy,.6*dx,.6*dy,width=.07,head_width=.3,head_length=.25,length_includes_head=True,color='#333',lw=0))
    un=[p for p,l in g['labs'].items() if l=='z3unhappy']
    a.plot([p[0] for p in un],[p[1] for p in un],'o',ms=2.6,color='#c0392b',mec='none',zorder=5)
    a.set_xlim(x0-.6,x1+.6); a.set_ylim(y0-.6,y1+.6); a.set_aspect('equal'); a.axis('off')
def arctic(key):
    d=D[key]; cells=sorted(d['cells']); n=d['n']; sh0=(n+1)%2 if key=='I' else n%2
    Xc=np.array([[c[0]+.5,c[1]+.5,1] for c in cells]); Yc=np.array([pcen(c) for c in cells]); Aff,*_=np.linalg.lstsq(Xc,Yc,rcond=None)
    th=np.linspace(0,2*np.pi,400); return np.c_[sh0+n/math.sqrt(2)*np.cos(th),n/math.sqrt(2)*np.sin(th),np.ones_like(th)]@Aff
EN={'II':F[5]['free']['E'],'forced':F[5]['forced']['E'],'I':30.0}
fig=plt.figure(figsize=(TWO,6.9)); TITLES=[]; HEAD=[]
gs=fig.add_gridspec(3,4,width_ratios=[1,1,1,0.42],wspace=0.1,hspace=0.36,left=0.0,right=0.995,top=0.935,bottom=0.012)
rows=(('II','edges along block boundaries, edge plaquettes free (ground state): Aztec II'),('forced','edges along block boundaries, every plaquette forced into a tetromino (comparison)'),('I','edges through blocks: Aztec I'))
for r,(key,lab) in enumerate(rows):
    L=['abc','def','ghi'][r]
    s=fig.add_subplot(gs[r,0]); gs_panel(s,key); TITLES.append((s,f'({L[0]}) vertex model, $12\\times12$, $E={G[(5,key)]["E"] if "E" in G[(5,key)] else EN[key]:.0f}J$'))
    a=fig.add_subplot(gs[r,1]); draw(a,key,'amp')
    if key!='forced': c=arctic(key); a.plot(c[:,0],c[:,1],color='white',lw=1.0,ls='--')
    TITLES.append((a,f'({L[1]}) $|\\Psi|/\\Psi_{{\\rm sat}}$, $48\\times48$'))
    b=fig.add_subplot(gs[r,2]); draw(b,key,'phase'); TITLES.append((b,f'({L[2]}) $\\arg\\Psi$, $48\\times48$'))
    HEAD.append((s,lab))
cax=fig.add_axes([0.905,0.70,0.014,0.17]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=amp_norm,cmap='viridis'),cax=cax,ticks=[0,.5,1])
cb.ax.tick_params(labelsize=6.5,length=2,width=.5); cb.outline.set_linewidth(.5); cb.set_label('$|\\Psi|/\\Psi_{\\rm sat}$',fontsize=7)
w=fig.add_axes([0.875,0.10,0.11,0.14],projection='polar'); wheel(w,fs=6); w.set_title('$\\arg\\Psi$',pad=7,fontsize=7.5)
fig.canvas.draw()
for ax_,t in TITLES:
    bb=ax_.get_position(); fig.text(bb.x0,bb.y1+0.008,t,ha='left',va='bottom',fontsize=7.5)
for ax_,t in HEAD:
    bb=ax_.get_position(); fig.text(0.005,bb.y1+0.032,t,fontsize=8,fontweight='bold',ha='left',va='bottom')
save(fig,'fig_axis')
