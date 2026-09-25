from pubstyle import *
import pickle, math, cmath, io, contextlib, sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pcolor import colour, wheel, outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, pcen
pl=sys.argv[1] if len(sys.argv)>1 else 'I'
M=pickle.load(open(f'exact/psi_{pl}.pkl','rb'))
def frame(D):
    cx,cy=D['centre']; Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]
    return sum(p[0] for p in Cc)/4, sum(p[1] for p in Cc)/4
def panel(a,D,P,mono):
    X,Y=frame(D); Wd=9.5; polys=[]; fcs=[]
    for c,z in P.items():
        px,py=pcen(c)
        if abs(px-X)>Wd+1 or abs(py-Y)>Wd+1: continue
        polys.append(Polygon(outline(list(rect(c))),closed=True)); fcs.append(colour(z))
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors='#999',linewidths=.15))
    for c,z in P.items():
        if (c[0]+c[1])%2: continue
        px,py=pcen(c)
        if abs(px-X)>Wd-.7 or abs(py-Y)>Wd-.7 or abs(z)<.03: continue
        u=z/abs(z)*0.8; a.arrow(px-u.real/2,py-u.imag/2,u.real,u.imag,head_width=.4,head_length=.35,length_includes_head=True,color='k',lw=.35)
    if mono:
        A,B=D['AB']
        for c in (A,B): a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#222',ec='k',lw=.5)); a.plot(*pcen(c),'*',color='w',ms=4.5)
    a.set_xlim(X-Wd,X+Wd); a.set_ylim(Y-Wd,Y+Wd); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    for s_ in a.spines.values(): s_.set_linewidth(.5)
fig=plt.figure(figsize=(TWO,3.55))
gs=fig.add_gridspec(2,5,width_ratios=[1,1,1,1,0.5],wspace=0.05,hspace=0.14,left=0.03,right=0.995,top=0.94,bottom=0.01)
for j,n in enumerate(sorted(M)):
    D=M[n]
    a=fig.add_subplot(gs[0,j]); panel(a,D,D['none'][0],False); a.set_title(f'({"abcd"[j]}) $n={n}$ ($n\\equiv{n%4}$)',loc='left',pad=2)
    b=fig.add_subplot(gs[1,j]); panel(b,D,D['mono'][0],True); b.set_title(f'({"efgh"[j]}) monomer pair',loc='left',pad=2)
    if j==0: a.set_ylabel('no defect',fontsize=8); b.set_ylabel('half vortices',fontsize=8)
w=fig.add_subplot(gs[0,4],projection='polar'); wheel(w,fs=6); w.set_title('$\\arg\\Psi$',pad=7,fontsize=7.5)
t=fig.add_subplot(gs[1,4]); t.axis('off'); t.text(0.05,0.8,'brightness:\n$|\\Psi|$\n(0 to 0.45)\n\narrows:\ndirection\nof $\\Psi$',fontsize=6.5,va='top',transform=t.transAxes)
save(fig,f'fig8_psi_{pl}')
