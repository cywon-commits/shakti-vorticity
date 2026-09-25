import pickle, math, cmath, numpy as np, glob
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect
L=48
def load(H,k): return pickle.load(open(f'hole/s_{L}_{H}_k{k}.pkl','rb'))
def ring_fourier(V,cell,d0,d1):
    a=0j; b=0j; n=0
    for i in range(L):
        for j in range(L):
            if not cell[i,j]: continue
            d=max(abs(i+.5-L/2),abs(j+.5-L/2))
            if d0<=d<d1:
                n+=1; a+=V[i,j]*cmath.exp(-1j*math.pi/2*(i+j)); b+=V[i,j]*cmath.exp(-1j*math.pi/2*(i-j))
    return a/n,b/n
base=load(0,0)
runs=[(8,0),(8,1),(8,-1),(10,0),(10,1)]
prof={}
for H,k in runs:
    try: D=load(H,k)
    except FileNotFoundError: continue
    rows=[]
    for d0 in np.arange(H/2,L/2,2):
        ah,bh=ring_fourier(D['vmean'],D['cell'],d0,d0+2); ab,bb=ring_fourier(base['vmean'],D['cell'],d0,d0+2)
        Ah=math.sqrt(2*(abs(ah)**2+abs(bh)**2)); Ab=math.sqrt(2*(abs(ab)**2+abs(bb)**2))
        align=((ah*ab.conjugate()+bh*bb.conjugate()).real)/(math.sqrt(abs(ah)**2+abs(bh)**2)*math.sqrt(abs(ab)**2+abs(bb)**2)+1e-12)
        rows.append((d0+1,Ah,Ab,align,abs(ah),abs(bh)))
    prof[(H,k)]=rows
    print(f'H={H} k={k}  Σv={sorted(set(D["tot"].tolist()))}')
    for r in rows: print('   ring d=%4.1f  A_hole=%.3f  A_nohole=%.3f  alignment cos=%+.2f   |alpha|=%.3f |beta|=%.3f'%r)
pickle.dump(prof,open('hole/prof.pkl','wb'))
# ---- figure: profiles + maps
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.6,.6)
fig=plt.figure(figsize=(20,11)); gs=fig.add_gridspec(2,4)
axp=fig.add_subplot(gs[0,0]); axa=fig.add_subplot(gs[1,0])
cols={(8,0):'#1b6ca8',(8,1):'#b0413e',(8,-1):'#e0892b',(10,0):'#3a9d5d',(10,1):'#8e44ad'}
base_rows=None
for key_,rows in prof.items():
    lab=f'hole {key_[0]}, sector k={key_[1]:+d}'
    axp.plot([r[0] for r in rows],[r[1] for r in rows],'o-',color=cols[key_],ms=4,label=lab)
    axa.plot([r[0] for r in rows],[r[3] for r in rows],'o-',color=cols[key_],ms=4,label=lab)
    base_rows=rows if key_==(8,0) else base_rows
R0=[(d0+1,*ring_fourier(base['vmean'],base['cell'],d0,d0+2)) for d0 in np.arange(0,L/2,2)]
axp.plot([r[0] for r in R0],[math.sqrt(2*(abs(r[1])**2+abs(r[2])**2)) for r in R0],'k--',label='no hole')
axp.set_xlabel('ring distance from centre (cells)'); axp.set_ylabel('local superstructure amplitude'); axp.legend(frameon=False,fontsize=8); axp.grid(alpha=.25)
axa.axhline(0,color='#999',lw=.8); axa.set_ylim(-1.05,1.05)
axa.set_xlabel('ring distance from centre (cells)'); axa.set_ylabel('alignment with the no-hole pattern (cos)'); axa.grid(alpha=.25)
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
panels=[((0,0),'no hole'),((8,0),'hole 8, k=0 (99.1%)'),((8,1),'hole 8, k=+1 (forced, 0.4%)'),((10,0),'hole 10, k=0 (50%)'),((10,1),'hole 10, k=+1 (50%)'),((8,-1),'hole 8, k=−1 (forced, 0.4%)')]
pos=[gs[0,1],gs[0,2],gs[0,3],gs[1,1],gs[1,2],gs[1,3]]
for (key_,title),g in zip(panels,pos):
    try: D=load(*key_)
    except FileNotFoundError: continue
    a=fig.add_subplot(g); V=D['vmean']; W=13
    for i in range(L):
        for j in range(L):
            if D['cell'][i,j] and abs(i+.5-L/2)<=W and abs(j+.5-L/2)<=W:
                a.add_patch(Polygon(outline(list(rect((i,j)))),closed=True,fc=cmapv(nv(V[i,j])),ec='none'))
    xs=[p[0] for i in range(L) for j in range(L) if abs(i+.5-L/2)<=W and abs(j+.5-L/2)<=W for p in rect((i,j))]
    ys=[p[1] for i in range(L) for j in range(L) if abs(i+.5-L/2)<=W and abs(j+.5-L/2)<=W for p in rect((i,j))]
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off'); a.set_title(title+' — Shakti view',fontsize=11)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=fig.axes[2:],shrink=.5,pad=.01); cb.set_label('⟨v⟩')
fig.suptitle('Square board L = 48 with a central hole: fixed winding sector = fixed rim vorticity',fontsize=13)
plt.savefig('outputs/hole_interference.png',dpi=140); plt.savefig('outputs/hole_interference.pdf')
