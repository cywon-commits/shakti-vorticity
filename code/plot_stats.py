import pickle, glob, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import *
PUR,NEU,YEL='#6a3d9a','#e6e6e6','#e8a33a'
cmap=LinearSegmentedColormap.from_list('py',[PUR,NEU,YEL]); norm=Normalize(-1,1)
seq=LinearSegmentedColormap.from_list('g',['#f4f4f4','#2b6a8f'])
def load(n,mode): return pickle.load(open(f'stats_{n}_{mode}.pkl','rb'))
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
def chi(c): return 1 if black(c) else -1
def shakti_map(ax,D,val,cm,nm,ring,title):
    for c in D['cells']:
        ax.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cm(nm(val(c))),ec='#999',lw=.4))
    for c in ring:
        r=rect(c); x=sum(p[0] for p in r)/2+.5; y=sum(p[1] for p in r)/2+.5
        ax.plot(x,y,'o',ms=5,mfc='none',mec='k',mew=1.2)
    xs=[p[0] for c in D['cells'] for p in rect(c)]; ys=[p[1] for c in D['cells'] for p in rect(c)]
    ax.set_xlim(min(xs)-.5,max(xs)+1.5); ax.set_ylim(min(ys)-.5,max(ys)+1.5)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=11)
def chess_map(ax,D,val,cm,nm,title,arctic=True):
    n=D['n']; cells=D['cells']
    xs=sorted({c[0] for c in cells}); ys=sorted({c[1] for c in cells})
    img=np.full((len(ys),len(xs),4),1.0)
    for c in cells: img[c[1]-ys[0],c[0]-xs[0]]=cm(nm(val(c)))
    ax.imshow(img,origin='lower',extent=(xs[0],xs[-1]+1,ys[0],ys[-1]+1),interpolation='nearest')
    if arctic:
        cx=(xs[0]+xs[-1]+1)/2; cy=(ys[0]+ys[-1]+1)/2; t=np.linspace(0,2*np.pi,300)
        ax.plot(cx+n/math.sqrt(2)*np.cos(t),cy+n/math.sqrt(2)*np.sin(t),color='k',lw=.8,ls='--')
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=11)
def ringof(D):
    n=D['n']; cells=D['cells']; sh=((n+1)%2,0)
    d=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n)]
    return analyse(d,cells)[2]
# ---------------- Figure 1: n=5 exact ----------------
D=load(5,'exact'); A,N=D['A'],D['N']; ring=ringof(D)
fig,ax=plt.subplots(1,3,figsize=(16,6))
shakti_map(ax[0],D,lambda c:A['sig'][c]/N,cmap,norm,ring,'⟨σ⟩  (colour of the short arm; free component → 0)')
shakti_map(ax[1],D,lambda c:A['tau'][c]/N,cmap,norm,ring,'⟨τ⟩ = ⟨σ·χ⟩  (χ = +1 horizontal, −1 vertical)')
shakti_map(ax[2],D,lambda c:A['conn'][c]/N,seq,Normalize(0,1),ring,'P(short arm connected to the ring)')
for a,cm,nm,lab in [(ax[0],cmap,norm,'⟨σ⟩'),(ax[1],cmap,norm,'⟨τ⟩'),(ax[2],seq,Normalize(0,1),'probability')]:
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=nm,cmap=cm),ax=a,shrink=.7,orientation='horizontal',pad=.02); cb.set_label(lab)
fig.suptitle('n = 5, exact average over all 32 768 tilings (placement I; ring spins fixed, ○ = frozen ring)',fontsize=12)
plt.tight_layout(); plt.savefig('outputs/ring_stats_n5_exact.png',dpi=200); plt.savefig('outputs/ring_stats_n5_exact.pdf')
# ---------------- Figure 2: large n, chessboard view ----------------
fig,ax=plt.subplots(2,3,figsize=(16,11))
for a,(n,m) in zip(ax[0],[(20,'3000'),(21,'3000'),(18,'1500')]):
    D=load(n,m); A,N=D['A'],D['N']
    chess_map(a,D,lambda c:A['tau'][c]/N,cmap,norm,f'n = {n} (n mod 4 = {n%4}):  ⟨τ⟩,  {N} samples')
for a,(n,m) in zip(ax[1],[(40,'1000'),(41,'1000')]):
    D=load(n,m); A,N=D['A'],D['N']
    chess_map(a,D,lambda c:A['tau'][c]/N,cmap,norm,f'n = {n} (n mod 4 = {n%4}):  ⟨τ⟩,  {N} samples')
D=load(40,'1000'); A,N=D['A'],D['N']
chess_map(ax[1][2],D,lambda c:A['conn'][c]/N,seq,Normalize(0,1),'n = 40:  P(short arm connected to the ring)')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax[:,:2].ravel().tolist(),shrink=.5,orientation='horizontal',pad=.03); cb.set_label('⟨τ⟩  (staggered colour; purple −1, yellow +1)')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=seq),ax=[ax[1][2]],shrink=.8,orientation='horizontal',pad=.03); cb.set_label('probability')
fig.suptitle('Chessboard view (dashed: arctic circle).  τ = σ·χ is constant along every connected chain, so it shows which ring segment a region is slaved to',fontsize=12)
plt.savefig('outputs/ring_stats_large_n.png',dpi=170); plt.savefig('outputs/ring_stats_large_n.pdf')
# ---------------- Figure 3: trends ----------------
rows=[]
for f in glob.glob('stats_*.pkl'):
    D=pickle.load(open(f,'rb'))
    if D['N']<500: continue
    g=D['glob']; fr=np.array([x['nconn']/x['nshort'] for x in g]); Mt=np.array([x['Mtau'] for x in g])
    rows.append((D['n'],fr.mean(),fr.std()/math.sqrt(len(g)),Mt.mean(),Mt.std()/math.sqrt(len(g)),np.mean([x['nshort']-x['nring'] for x in g])))
rows.sort()
nn=np.array([r[0] for r in rows])
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
ax[0].errorbar(nn,[r[1] for r in rows],yerr=[r[2] for r in rows],fmt='o-',color='#2b6a8f',ms=6,lw=2,capsize=2)
ax[0].set_xlabel('n'); ax[0].set_ylabel('fraction of short arms connected to the ring'); ax[0].set_ylim(0,1); ax[0].grid(alpha=.25)
for k,mk,colr in [(2,'s','#c0392b'),(0,'o','#555'),(1,'^','#555'),(3,'v','#555')]:
    sel=[r for r in rows if r[0]%4==k]
    ax[1].errorbar([r[0] for r in sel],[abs(r[3]) for r in sel],yerr=[r[4] for r in sel],fmt=mk+('-' if k==2 else ''),color=colr,ms=7,lw=2,capsize=2,
                   mfc=colr if k==2 else 'white',label=f'n ≡ {k} (mod 4)')
ax[1].set_xlabel('n'); ax[1].set_ylabel('|⟨M_τ⟩| = |⟨Σ_interior τ⟩|'); ax[1].grid(alpha=.25); ax[1].legend(frameon=False)
plt.tight_layout(); plt.savefig('outputs/ring_stats_trends.png',dpi=200); plt.savefig('outputs/ring_stats_trends.pdf')
for r in rows: print(r)
