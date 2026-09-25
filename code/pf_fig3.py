from pubstyle import *
import random, pickle
from matplotlib.patches import Polygon
src=open('plot_vort.py').read().split("fig,ax=plt.subplots(2,3")[0]
src=src.replace("import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt","")
exec(src)
import matplotlib.pyplot as plt
nv=Normalize(-1,1)
def draw_ex(ax,n,seed):
    cells,sh=aztec(n); rng=random.Random(seed)
    doms=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
    for k,d in enumerate(doms):
        A=info2(d); s=v[k]*ref
        ax.add_patch(Polygon(outline(list(A['S'])),closed=True,fc=OR if s>0 else BL,ec='white',lw=.3))
        ax.add_patch(Polygon(outline(list(A['L'])),closed=True,fc=ORl if s>0 else BLl,ec='white',lw=.3))
        ax.add_patch(Polygon(outline(list(A['S'])+list(A['L'])),closed=True,fill=False,ec='k',lw=.7))
        (p,q)=A['defect']; ax.plot([p[0],q[0]],[p[1],q[1]],color='#c00000',lw=1.8,solid_capstyle='butt',zorder=6)
    xs=[p[0] for c in cells for p in rect(c)]; ys=[p[1] for c in cells for p in rect(c)]
    ax.set_xlim(min(xs)-.3,max(xs)+1.3); ax.set_ylim(min(ys)-.3,max(ys)+1.3); ax.set_aspect('equal'); ax.axis('off')
    return sum(v[k]*ref for k in range(len(doms)))
def avg(ax,D):
    A,N=D['A'],D['N']
    for c in D['cells']: ax.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=CMAPV(nv(A['v'][c]/N)),ec='#999',lw=.25))
    xs=[p[0] for c in D['cells'] for p in rect(c)]; ys=[p[1] for c in D['cells'] for p in rect(c)]
    ax.set_xlim(min(xs)-.3,max(xs)+1.3); ax.set_ylim(min(ys)-.3,max(ys)+1.3); ax.set_aspect('equal'); ax.axis('off')
fig,ax=plt.subplots(1,4,figsize=(TWO,1.95),gridspec_kw=dict(wspace=0.04,left=0,right=0.935,top=0.9,bottom=0))
V5=draw_ex(ax[0],5,3); V6=draw_ex(ax[1],6,3)
ax[0].set_title(f'(a) $n=5$, one tiling ($V={V5:+d}$)',pad=1); ax[1].set_title(f'(b) $n=6$, one tiling ($V={V6:+d}$)',pad=1)
avg(ax[2],load(5,'exact')); avg(ax[3],load(6,'20000'))
ax[2].set_title('(c) $n=5$: $\\langle v\\rangle$ (exact)',pad=1); ax[3].set_title('(d) $n=6$: $\\langle v\\rangle$ (MC)',pad=1)
small_cbar(fig,plt.cm.ScalarMappable(norm=nv,cmap=CMAPV),ax[3],label='$\\langle v\\rangle$',size=('5%','70%'),pad=0.03,ticks=[-1,0,1])
save(fig,'fig3_examples')
