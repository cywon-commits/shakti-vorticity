import pickle, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect, pcen
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
rows=pickle.load(open('sq/rows.pkl','rb'))
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
# ---- scaling
fig,ax=plt.subplots(1,3,figsize=(17,4.8))
st={('black',0):('#1b6ca8','o','full'),('black',2):('#1b6ca8','o','none'),('white',0):('#e0892b','s','full'),('white',2):('#e0892b','s','none')}
for (cls,par),(c,m,fill) in st.items():
    R=[r for r in rows if r['cls']==cls and r['L']%4==par]
    lab=f'corner {cls}, L ≡ {par} (mod 4)'
    kw=dict(color=c,ms=7,mew=1.5,mfc=c if fill=='full' else 'white',ls='',marker=m,label=lab)
    ax[0].plot([r['L'] for r in R],[r['amp'] for r in R],**kw)
    ax[1].plot([r['L'] for r in R],[r['var'] for r in R],**kw)
    ax[2].plot([r['L'] for r in R],[abs(r['tot'][0]) for r in R],**kw)
LL=np.linspace(8,48,50); ax[0].plot(LL,0.42*(LL/16)**-0.25,'k--',lw=1.2,label='∝ L^(−1/4)')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xticks([8,16,32,48]); ax[0].set_xticklabels(['8','16','32','48'])
ax[0].set_yticks([.2,.3,.4,.5,.6]); ax[0].set_yticklabels(['0.2','0.3','0.4','0.5','0.6'])
ax[0].set_xlabel('L'); ax[0].set_ylabel('2×2-superblock amplitude of ⟨v⟩ (central disk)'); ax[0].legend(frameon=False,fontsize=8)
ax[1].set_xscale('log'); ax[1].set_xticks([8,16,32,48]); ax[1].set_xticklabels(['8','16','32','48']); ax[1].set_xlabel('L'); ax[1].set_ylabel('Var h(centre)')
ax[2].set_xlabel('L'); ax[2].set_ylabel('|Σ v|  (same for every tiling)'); ax[2].set_yticks([0,2])
for a in ax: a.grid(alpha=.25)
plt.tight_layout(); plt.savefig('outputs/square_scaling.png',dpi=190); plt.savefig('outputs/square_scaling.pdf')
# ---- maps: chessboard and Shakti, L = 32 and 34, both classes
fig,ax=plt.subplots(2,4,figsize=(22,11.5))
for j,(L,sh) in enumerate([(32,'00'),(32,'10'),(34,'00'),(34,'10')]):
    D=pickle.load(open(f'sq/sq_{L}_{sh}.pkl','rb')); V=D['vmean']; sx,sy=D['shift']
    cls='black' if (sx+sy)%2==0 else 'white'
    ax[0][j].imshow(V.T,origin='lower',cmap=cmapv,norm=nv,extent=(sx,sx+L,sy,sy+L))
    ax[0][j].set_title(f'L = {L}, corner {cls}: chessboard ⟨v⟩',fontsize=12); ax[0][j].set_xticks([]); ax[0][j].set_yticks([])
    a=ax[1][j]; xs=[];ys=[]
    for i in range(L):
        for k in range(L):
            c=(sx+i,sy+k); r=list(rect(c)); xs+= [p[0] for p in r]; ys+=[p[1] for p in r]
            a.add_patch(Polygon(outline(r),closed=True,fc=cmapv(nv(V[i,k])),ec='none'))
    a.set_xlim(min(xs),max(xs)+1); a.set_ylim(min(ys),max(ys)+1); a.set_aspect('equal'); a.axis('off')
    a.set_title(f'L = {L}, corner {cls}: Shakti ⟨v⟩ per plaquette',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.ravel().tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩')
fig.suptitle('Square L×L board (uniform tilings, Monte Carlo): whole board, chessboard view (top) and Shakti view (bottom)',fontsize=13)
plt.savefig('outputs/square_maps.png',dpi=140); plt.savefig('outputs/square_maps.pdf')
# ---- Shakti centre zoom (like the Aztec figure)
fig,ax=plt.subplots(2,2,figsize=(13,13)); ax=ax.ravel()
for a,(L,sh) in zip(ax,[(32,'00'),(32,'10'),(34,'00'),(34,'10')]):
    D=pickle.load(open(f'sq/sq_{L}_{sh}.pkl','rb')); V=D['vmean']; sx,sy=D['shift']; cx,cy=sx+L//2,sy+L//2
    C=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in C)/4; Y=sum(p[1] for p in C)/4; W=9
    for i in range(L):
        for k in range(L):
            c=(sx+i,sy+k); px,py=pcen(c)
            if abs(px-X)<=W and abs(py-Y)<=W: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cmapv(nv(V[i,k])),ec='#555',lw=.5))
    a.plot(X,Y,'k+',ms=14,mew=2); a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'L = {L}, corner {"black" if (sx+sy)%2==0 else "white"}: ⟨v⟩ near the centre',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩')
fig.suptitle('Square board, centre on the Shakti lattice (+ = centre)',fontsize=13)
plt.savefig('outputs/square_centre_shakti.png',dpi=160); plt.savefig('outputs/square_centre_shakti.pdf')
