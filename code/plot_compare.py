import pickle, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect, pcen
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
R1=[r for r in pickle.load(open('fast/rows.pkl','rb')) if r['N']>=600]
R2=[r for r in pickle.load(open('fast2/rows.pkl','rb')) if r['N']>=600]
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
for R,lab,c,m in ((R1,'placement I','#1b6ca8','o'),(R2,'placement II','#e0892b','s')):
    ax[0].plot([r['n'] for r in R],[r['mv'] for r in R],m,color=c,mfc='none' if lab.endswith('II') else c,ms=7,mew=1.5,label=lab)
    ax[1].plot([r['n'] for r in R],[abs(r['tot']) for r in R],m,color=c,mfc='none' if lab.endswith('II') else c,ms=7,mew=1.5,label=lab)
nn=np.linspace(6,80,100); ax[0].plot(nn,0.35*(nn/20)**-0.25,'k--',lw=1.2,label='∝ n^(−1/4)')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xticks([6,10,20,40,80]); ax[0].set_xticklabels(['6','10','20','40','80'])
ax[0].set_yticks([.2,.3,.4,.5]); ax[0].set_yticklabels(['0.2','0.3','0.4','0.5'])
ax[0].set_xlabel('n'); ax[0].set_ylabel('2×2-superblock amplitude of ⟨v⟩ (centre)'); ax[0].legend(frameon=False)
ax[1].set_xlabel('n'); ax[1].set_ylabel('|Σ v|  (same for every tiling)'); ax[1].set_yticks([0,2]); ax[1].legend(frameon=False)
for a in ax: a.grid(alpha=.25)
plt.tight_layout(); plt.savefig('outputs/placement_I_vs_II_scaling.png',dpi=190); plt.savefig('outputs/placement_I_vs_II_scaling.pdf')
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
# chessboard + Shakti centre maps, placement II
fig,ax=plt.subplots(2,4,figsize=(21,10.5))
for j,n in enumerate([36,37,38,39]):
    D=pickle.load(open(f'fast2/f_{n}_1500.pkl','rb')); x0,y0=D['x0'],D['y0']; cx,cy=D['center']; V=D['vmean']; w=8
    img=np.full((2*w,2*w),np.nan)
    for i in range(-w,w):
        for k in range(-w,w):
            c=(cx+i,cy+k)
            if c in D['cells']: img[k+w,i+w]=V[c[0]-x0,c[1]-y0]
    ax[0][j].imshow(img,origin='lower',cmap=cmapv,norm=nv,extent=(cx-w,cx+w,cy-w,cy+w)); ax[0][j].plot(cx,cy,'k+',ms=12,mew=2)
    ax[0][j].set_title(f'II, n = {n} (≡{n%4}): chessboard ⟨v⟩',fontsize=12); ax[0][j].set_xticks([]); ax[0][j].set_yticks([])
    C=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in C)/4; Y=sum(p[1] for p in C)/4; W=9
    a=ax[1][j]
    for c in D['cells']:
        px,py=pcen(c)
        if abs(px-X)<=W and abs(py-Y)<=W: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cmapv(nv(V[c[0]-x0,c[1]-y0])),ec='#555',lw=.5))
    a.plot(X,Y,'k+',ms=12,mew=2); a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'II, n = {n}: Shakti ⟨v⟩ per plaquette',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.ravel().tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩')
fig.suptitle('Placement II, centre of the Aztec diamond (+ = centre; 1500 tilings per n)',fontsize=13)
plt.savefig('outputs/placement_II_centre.png',dpi=150); plt.savefig('outputs/placement_II_centre.pdf')
