import pickle, math, cmath, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from core import rect
from plot_psi import colour, wheel
from plot_defect_maps import outline
D=pickle.load(open('exact/phase_map.pkl','rb'))
names=['Aztec I n=24','square L=40','rectangle 56x28']
fig=plt.figure(figsize=(20,13)); gs=fig.add_gridspec(3,4,width_ratios=[1,1,1,.6])
for j,name in enumerate(names):
    d=D[name]; P,hc=d['P'],d['hc']
    for i,(title,f) in enumerate([('exact arg Ψ (Kasteleyn)',lambda c:colour(P[c]/max(abs(P[c]),1e-9)*0.45)),
                                  ('prediction π⟨h⟩/4 (exact mean height)',lambda c:colour(0.45*cmath.exp(1j*math.pi*hc[c]/4))),
                                  ('|arg Ψ − π⟨h⟩/4|',None)]):
        a=fig.add_subplot(gs[i,j]); xs=[];ys=[]
        for c in d['cells']:
            if c not in P: continue
            r=list(rect(c)); xs+=[p[0] for p in r]; ys+=[p[1] for p in r]
            if f is None:
                dv=abs(((cmath.phase(P[c])-math.pi*hc[c]/4+math.pi)%(2*math.pi))-math.pi)
                fc=plt.cm.Greys(min(dv/0.4,1))
            else: fc=f(c)
            a.add_patch(Polygon(outline(r),closed=True,fc=fc,ec='none'))
        a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off')
        a.set_title(f'{name}: {title}',fontsize=10.5)
wheel(fig.add_subplot(gs[0,3],projection='polar'))
cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,.4),cmap='Greys'),ax=fig.add_subplot(gs[2,3]),fraction=.5); cb.set_label('phase error (rad)'); fig.axes[-2].axis('off')
t=fig.add_subplot(gs[1,3]); t.axis('off')
t.text(0,1,"Phase map: arg Ψ(x) = π⟨h(x)⟩/4\n\n• square / rectangle: rms error 0.006 rad\n• Aztec: exact at the centre (symmetry),\n  < 0.1 rad in the inner half of the disk,\n  larger only near the arctic circle\n\nCentre, 90°-symmetric region:\n⟨h(c)⟩ = [h(p) + h(Rp)]/2 exactly\n(p any boundary vertex, R = 90° rotation\nabout c, which swaps colours: h → c₀ − h)\n\nAztec I:  ⟨h(c)⟩ = 3 − n (n even), 1 − n (odd)\nAztec II: ⟨h(c)⟩ = n + 2 (n even), n (odd)\n(normalisation h = r(p) at the bottom corner)",va='top',fontsize=10.5)
fig.suptitle('Kenyon–Okounkov phase map: the vorticity phase is set by the mean height function',fontsize=14)
plt.tight_layout(rect=(0,0,1,.96)); plt.savefig('outputs/phase_map.png',dpi=120); plt.savefig('outputs/phase_map.pdf'); print('ok')
