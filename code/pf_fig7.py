from pubstyle import *
import pickle, math, cmath, io, contextlib, numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pcolor import colour, wheel, outline
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, pcen
D=pickle.load(open('exact/phase_map.pkl','rb'))
err_norm=Normalize(0,0.4); amp_norm=Normalize(0,1)
def draw(a,key,mode):
    d=D[key]; P,hc=d['P'],d['hc']; polys=[]; fcs=[]; xs=[];ys=[]
    for c in d['cells']:
        if c not in P: continue
        r=list(rect(c)); xs+=[p[0] for p in r]; ys+=[p[1] for p in r]; polys.append(Polygon(outline(r),closed=True))
        if mode=='phase': fcs.append(colour(P[c]/max(abs(P[c]),1e-9)*0.45))
        elif mode=='amp': fcs.append(plt.cm.viridis(amp_norm(abs(P[c])/PSI_SAT)))
        else: fcs.append(plt.cm.Greys(err_norm(abs(((cmath.phase(P[c])-math.pi*hc[c]/4+math.pi)%(2*math.pi))-math.pi))))
    a.add_collection(PatchCollection(polys,facecolors=fcs,edgecolors=fcs,linewidths=.3))
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off')
# --- Aztec n=24 (placement I): saturation value of |Psi| in the frozen corners, arctic circle in Shakti coordinates
n=24; AZ=D['Aztec I n=24']; sh0=(n+1)%2
def rad(c): return math.hypot(c[0]+.5-sh0,c[1]+.5)/(n/math.sqrt(2))     # distance from the centre in units of the arctic radius
PSI_SAT=np.median([abs(v) for c,v in AZ['P'].items() if rad(c)>1.15])     # fully ordered (frozen) value
cells=sorted(AZ['P']); X=np.array([[c[0]+.5,c[1]+.5,1] for c in cells]); Y=np.array([pcen(c) for c in cells])
Aff,*_=np.linalg.lstsq(X,Y,rcond=None)                                     # chessboard -> Shakti (affine)
th=np.linspace(0,2*np.pi,400); circ=np.c_[sh0+n/math.sqrt(2)*np.cos(th),n/math.sqrt(2)*np.sin(th),np.ones_like(th)]@Aff
fig=plt.figure(figsize=(TWO,3.75)); TITLES=[]
gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,1.0],height_ratios=[1,1],wspace=0.22,hspace=0.48,left=0.0,right=0.99,top=0.93,bottom=0.03)
a=fig.add_subplot(gs[0,0]); draw(a,'Aztec I n=24','phase'); TITLES.append((a,'(a) Aztec, $n=24$: $\\arg\\Psi$'))
b=fig.add_subplot(gs[0,1]); draw(b,'Aztec I n=24','amp'); TITLES.append((b,'(b) amplitude $|\\Psi|/\\Psi_{\\rm sat}$'))
b.plot(circ[:,0],circ[:,1],color='white',lw=1.0,ls='--')
cax=b.inset_axes([1.03,0.08,0.05,0.5]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=amp_norm,cmap='viridis'),cax=cax,ticks=[0,.5,1])
cb.ax.tick_params(labelsize=6,length=2,width=.5,pad=1); cb.outline.set_linewidth(.5)
# radial profile
c_=fig.add_subplot(gs[0,2:]); rr=np.array([rad(c) for c in cells]); aa=np.array([abs(AZ['P'][c])/PSI_SAT for c in cells])
c_.plot(rr,aa,'.',ms=1.6,color=CAT[0],alpha=.5,rasterized=True)
edges=np.arange(0,1.42,0.05); mid=(edges[1:]+edges[:-1])/2; mean=[aa[(rr>=lo)&(rr<hi)].mean() if ((rr>=lo)&(rr<hi)).any() else np.nan for lo,hi in zip(edges,edges[1:])]
c_.plot(mid,mean,'-',color='k',lw=1.1,label='ring average')
c_.axvline(1,color='#888',ls='--',lw=.8); c_.text(1.02,0.15,'arctic\ncircle',fontsize=6.5,color='#555',va='bottom')
c_.text(0.03,0.93,'liquid: algebraic order',fontsize=7,transform=c_.transAxes,va='top'); c_.text(0.97,0.82,'frozen',fontsize=7,transform=c_.transAxes,va='top',ha='right')
c_.set_xlim(0,1.42); c_.set_ylim(0,1.08); c_.set_xlabel('distance from centre / arctic radius $n/\\sqrt{2}$'); c_.set_ylabel('$|\\Psi|/\\Psi_{\\rm sat}$')
c_.legend(loc='center left',fontsize=6.5); TITLES.append((c_,'(c) radial profile'))
fig.canvas.draw(); pc=c_.get_position(); c_.set_position([pc.x0+0.075,pc.y0,pc.width-0.075,pc.height])
d=fig.add_subplot(gs[1,0]); draw(d,'Aztec I n=24','err'); TITLES.append((d,'(d) $|\\arg\\Psi-\\pi\\langle h\\rangle/4|$'))
cax=d.inset_axes([1.03,0.08,0.05,0.5]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=err_norm,cmap='Greys'),cax=cax,ticks=[0,.2,.4])
cb.ax.tick_params(labelsize=6,length=2,width=.5,pad=1); cb.outline.set_linewidth(.5); cb.set_label('rad',fontsize=6.5,labelpad=1)
e=fig.add_subplot(gs[1,1]); draw(e,'square L=40','phase'); TITLES.append((e,'(e) square, $L=40$'))
f=fig.add_subplot(gs[1,2]); draw(f,'rectangle 56x28','phase'); TITLES.append((f,'(f) rectangle, $56\\times28$'))
w=fig.add_subplot(gs[1,3],projection='polar'); wheel(w,fs=6); w.set_title('$\\arg\\Psi$',pad=7,fontsize=7.5)
pos=w.get_position(); w.set_position([pos.x0+pos.width*0.2,pos.y0+pos.height*0.15,pos.width*0.6,pos.height*0.6])
fig.canvas.draw()
for ax_,t in TITLES:
    bb=ax_.get_position(); fig.text(bb.x0,bb.y1+0.012,t,ha='left',va='bottom',fontsize=8)
print('Psi_sat = %.4f'%PSI_SAT)
save(fig,'fig7_phasemap')
