from pubstyle import *
import pickle, math, cmath, io, contextlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pcolor import colour, wheel, outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
D=pickle.load(open('exact/phase_map.pkl','rb'))
err_norm=Normalize(0,0.4)
def draw(a,key,mode):
    d=D[key]; P,hc=d['P'],d['hc']; polys=[]; fcs=[]; xs=[];ys=[]
    for c in d['cells']:
        if c not in P: continue
        r=list(rect(c)); xs+=[p[0] for p in r]; ys+=[p[1] for p in r]; polys.append(Polygon(outline(r),closed=True))
        if mode=='phase': fcs.append(colour(P[c]/max(abs(P[c]),1e-9)*0.45))
        else: fcs.append(plt.cm.Greys(err_norm(abs(((cmath.phase(P[c])-math.pi*hc[c]/4+math.pi)%(2*math.pi))-math.pi))))
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors=fcs,linewidths=.3))
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off')
fig=plt.figure(figsize=(TWO,1.95)); TITLES=[]
gs=fig.add_gridspec(1,5,width_ratios=[1,1,1,1.3,0.42],wspace=0.18,left=0,right=0.995,top=0.86,bottom=0.02)
a=fig.add_subplot(gs[0]); draw(a,'Aztec I n=24','phase'); a.set_anchor('C'); TITLES.append((a,'(a) Aztec, $n=24$: $\\arg\\Psi$'))
b=fig.add_subplot(gs[1]); draw(b,'Aztec I n=24','err'); b.set_anchor('C'); TITLES.append((b,'(b) $|\\arg\\Psi-\\pi\\langle h\\rangle/4|$'))
cax=b.inset_axes([1.03,0.08,0.05,0.42]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=err_norm,cmap='Greys'),cax=cax,ticks=[0,.2,.4])
cb.ax.tick_params(labelsize=6,length=2,width=.5,pad=1); cb.outline.set_linewidth(.5); cb.set_label('rad',fontsize=6.5,labelpad=1)
c=fig.add_subplot(gs[2]); draw(c,'square L=40','phase'); c.set_anchor('C'); TITLES.append((c,'(c) square, $L=40$'))
d=fig.add_subplot(gs[3]); draw(d,'rectangle 56x28','phase'); d.set_anchor('C'); TITLES.append((d,'(d) rectangle, $56\\times28$'))
w=fig.add_subplot(gs[4],projection='polar'); wheel(w,fs=6); w.set_title('$\\arg\\Psi$',pad=7,fontsize=7.5)
fig.canvas.draw()
for ax_,t in TITLES:
    bb=ax_.get_position(); fig.text(bb.x0,0.9,t,ha='left',va='bottom',fontsize=8)
save(fig,'fig7_phasemap')
