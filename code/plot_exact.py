import pickle, math, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
out=pickle.load(open('exact/scan_8_34.pkl','rb'))+pickle.load(open('exact/scan_35_58.pkl','rb'))
out=[r for r in out if r['lnerr']<1e-5]
fig,ax=plt.subplots(1,3,figsize=(18,5))
cols={0:'#1b6ca8',1:'#e0892b',2:'#3a9d5d',3:'#b0413e'}; mk={0:'o',1:'s',2:'^',3:'D'}
for k in range(4):
    R=[r for r in out if r['n']%4==k]
    ax[0].plot([r['n'] for r in R],[r['A'] for r in R],mk[k],color=cols[k],ms=6,label=f'n ≡ {k} (mod 4)')
    if k in (0,3): ax[1].plot([r['n'] for r in R],[abs(r['eps_c']) for r in R],mk[k],color=cols[k],ms=6,label=f'n ≡ {k}')
nn=np.linspace(8,60,100)
ax[0].plot(nn,0.35*(nn/20)**-0.25,'k--',lw=1,label='∝ n^(−1/4)')
m=np.linspace(2,15,100); ax[1].plot(4*m,0.7670*m**-0.25,'k--',lw=1,label='0.767·m^(−1/4), m = ⌈n/4⌉')
for a in ax[:2]:
    a.set_xscale('log'); a.set_yscale('log'); a.set_xticks([8,16,32,58]); a.set_xticklabels(['8','16','32','58'])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.grid(alpha=.25); a.legend(frameon=False,fontsize=9); a.set_xlabel('n')
ax[0].set_yticks([.28,.32,.36,.4]); ax[0].set_yticklabels(['0.28','0.32','0.36','0.40']); ax[0].set_ylabel('exact local amplitude (cells within r ≤ 4 of the centre)')
ax[1].set_yticks([.4,.5,.6,.7]); ax[1].set_yticklabels(['0.4','0.5','0.6','0.7']); ax[1].set_ylabel('exact |⟨ε(centre)⟩|  (0 exactly for n ≡ 1, 2)')
E={r['n']//4:abs(r['eps_c']) for r in out if r['n']%4==0}; ms=sorted(E)
ax[2].plot([mm for mm in ms[1:]],[math.log(E[mm-1]/E[mm])/math.log(mm/(mm-1)) for mm in ms[1:] if mm-1 in E] if all(mm-1 in E for mm in ms[1:]) else [],'o-',color='#1b6ca8')
ax[2].axhline(0.25,color='k',ls='--',lw=1); ax[2].set_xlabel('m = n/4'); ax[2].set_ylabel('local exponent  −d ln|⟨ε⟩| / d ln m'); ax[2].grid(alpha=.25); ax[2].set_ylim(0.235,0.255)
plt.tight_layout(); plt.savefig('outputs/exact_scaling.png',dpi=180); plt.savefig('outputs/exact_scaling.pdf')
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
M=pickle.load(open('exact/maps.pkl','rb'))
fig,ax=plt.subplots(2,4,figsize=(21,10.5))
for j,n in enumerate([36,37,38,39]):
    D=M[n]; V=D['V']; cx,cy=D['centre']; w=8
    img=np.full((2*w,2*w),np.nan)
    for i in range(-w,w):
        for k in range(-w,w):
            c=(cx+i,cy+k)
            if c in V: img[k+w,i+w]=V[c]
    ax[0][j].imshow(img,origin='lower',cmap=cmapv,norm=nv,extent=(cx-w,cx+w,cy-w,cy+w)); ax[0][j].plot(cx,cy,'k+',ms=12,mew=2)
    ax[0][j].set_title(f'n = {n} (n≡{n%4}): exact ⟨v⟩, chessboard',fontsize=12); ax[0][j].set_xticks([]); ax[0][j].set_yticks([])
    a=ax[1][j]; from core import pcen
    Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in Cc)/4; Y=sum(p[1] for p in Cc)/4; W=9
    for c,val in V.items():
        px,py=pcen(c)
        if abs(px-X)<=W and abs(py-Y)<=W: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cmapv(nv(val)),ec='#555',lw=.5))
    a.plot(X,Y,'k+',ms=12,mew=2); a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'n = {n}: exact ⟨v⟩, Shakti view',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.ravel().tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩ (exact)')
fig.suptitle('Exact vorticity at the centre of the Aztec diamond (placement I), Kasteleyn with a sign seam — no sampling noise',fontsize=13)
plt.savefig('outputs/exact_centre_maps.png',dpi=150); plt.savefig('outputs/exact_centre_maps.pdf')
# full diamond n=24 and 39
fig,ax=plt.subplots(1,2,figsize=(15,7.5))
for a,n in zip(ax,(24,39)):
    D=M[n]; V=D['V']; xs=[];ys=[]
    for c,val in V.items():
        r=list(rect(c)); xs+=[p[0] for p in r]; ys+=[p[1] for p in r]
        a.add_patch(Polygon(outline(r),closed=True,fc=cmapv(nv(val)),ec='none'))
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off'); a.set_title(f'n = {n}: exact ⟨v⟩ on the whole diamond (Shakti view)',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩ (exact)')
plt.savefig('outputs/exact_full_maps.png',dpi=150); plt.savefig('outputs/exact_full_maps.pdf')
