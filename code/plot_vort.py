import random, pickle, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon, FancyArrow
exec(open('vort_rule.py').read().split('doms=')[0])
from vort_tiling import signed_graph, frustration
_h=open('plot_stats.py').read().split('# ---------------- Figure 1')[0].split('from core import *')[1]
exec(_h.split('def ringof')[0])
def chi(c): return 1 if black(c) else -1
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-1,1)
def chi(c): return 1 if black(c) else -1
def load(n,m): return pickle.load(open(f'vstats_{n}_{m}.pkl','rb'))
def draw_example(ax,n,seed,title):
    cells,sh=aztec(n); rng=random.Random(seed)
    doms=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
    for k,d in enumerate(doms):
        A=info2(d); s=v[k]*ref
        base='#e8743a' if s>0 else '#2f6db5'; light='#f6c6a8' if s>0 else '#aac6ea'
        ax.add_patch(Polygon(outline(list(A['S'])),closed=True,fc=base,ec='white',lw=.6))
        ax.add_patch(Polygon(outline(list(A['L'])),closed=True,fc=light,ec='white',lw=.6))
        ax.add_patch(Polygon(outline(list(A['S'])+list(A['L'])),closed=True,fill=False,ec='k',lw=1.4))
        (p,q)=A['defect']; ax.plot([p[0],q[0]],[p[1],q[1]],color='#c00000',lw=3.2,solid_capstyle='butt',zorder=6)
        # circulation mark
        cx=sum(x[0]+.5 for x in A['S'])/2; cy=sum(x[1]+.5 for x in A['S'])/2
        ax.text(cx,cy,'↺' if s>0 else '↻',ha='center',va='center',fontsize=11,color='white',fontweight='bold')
    xs=[p[0] for c in cells for p in rect(c)]; ys=[p[1] for c in cells for p in rect(c)]
    ax.set_xlim(min(xs)-.5,max(xs)+1.5); ax.set_ylim(min(ys)-.5,max(ys)+1.5); ax.set_aspect('equal'); ax.axis('off')
    M=sum(v[k]*ref for k in range(len(doms)))
    ax.set_title(f'{title}\nΣ v = {M:+d}',fontsize=11)
def shakti_avg(ax,D,title):
    A,N=D['A'],D['N']
    for c in D['cells']:
        ax.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cmapv(nv(A['v'][c]/N)),ec='#aaa',lw=.4))
    xs=[p[0] for c in D['cells'] for p in rect(c)]; ys=[p[1] for c in D['cells'] for p in rect(c)]
    ax.set_xlim(min(xs)-.5,max(xs)+1.5); ax.set_ylim(min(ys)-.5,max(ys)+1.5); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=11)
fig,ax=plt.subplots(2,3,figsize=(18,12))
draw_example(ax[0][0],5,3,'n = 5: one tiling, vorticity from the island rule')
draw_example(ax[1][0],6,3,'n = 6: one tiling')
draw_example(ax[0][1],5,8,'n = 5: another tiling')
draw_example(ax[1][1],6,8,'n = 6: another tiling')
shakti_avg(ax[0][2],load(5,'exact'),'n = 5: ⟨v⟩ per plaquette (exact, 32 768 tilings)')
shakti_avg(ax[1][2],load(6,'20000'),'n = 6: ⟨v⟩ per plaquette (20 000 samples)')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax[:,2].tolist(),shrink=.5,orientation='horizontal',pad=.03); cb.set_label('⟨v⟩   (blue −1 = clockwise,  orange +1 = counter-clockwise)')
fig.text(.5,.015,'dark = strong (short) arm, light = weak (long) arm, red bar = the defect island of each tetromino (red arrow in the template).  Boundary ring vorticities are fixed.',ha='center',fontsize=11)
plt.savefig('outputs/vorticity_n5_n6.png',dpi=160); plt.savefig('outputs/vorticity_n5_n6.pdf')
# large n
fig,ax=plt.subplots(1,4,figsize=(22,6.5))
for a,(n,m) in zip(ax,[(20,'500'),(21,'500'),(30,'500'),(40,'500')]):
    D=load(n,m); A,N=D['A'],D['N']
    chess_map(a,D,lambda c:A['v'][c]/N,cmapv,nv,f'n = {n}: ⟨v⟩ ({N} samples, chessboard view)')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.tolist(),shrink=.5,orientation='horizontal',pad=.04); cb.set_label('⟨v⟩')
plt.savefig('outputs/vorticity_large_n.png',dpi=150); plt.savefig('outputs/vorticity_large_n.pdf')
