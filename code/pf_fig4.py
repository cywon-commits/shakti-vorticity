from pubstyle import *
import pickle, numpy as np, io, contextlib
from matplotlib.patches import Polygon
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, pcen
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
nv=Normalize(-.5,.5)
# exact <v> windows around the centre, from maps2.py (placements I and II)
MAPS={pl:pickle.load(open(f'exact/maps_{pl}.pkl','rb')) for pl in ('I','II')}
def centre_height(n,pl):
    """exact <h(c)> from 90-degree symmetry (Sec. V C)"""
    if pl=='I': return 3-n if n%2==0 else 1-n
    return n+2 if n%2==0 else n
def pattern(n,pl):
    h=centre_height(n,pl)
    if h%2: return 'pinwheel (H+V)'
    return 'facing (H only)' if h%4==2 else 'facing (V only)'
fig,ax=plt.subplots(3,4,figsize=(TWO,5.25),gridspec_kw=dict(wspace=0.05,hspace=0.17,left=0.035,right=0.93,top=0.965,bottom=0.0))
NS=[36,37,38,39]; lab='abcdefghijkl'
for j,n in enumerate(NS):
    D=MAPS['I'][n]; V=D['V']; cx,cy=D['centre']; w=8
    img=np.full((2*w,2*w),np.nan)
    for i in range(-w,w):
        for k in range(-w,w):
            if (cx+i,cy+k) in V: img[k+w,i+w]=V[(cx+i,cy+k)]
    a=ax[0][j]; a.imshow(img,origin='lower',cmap=CMAPV,norm=nv,extent=(cx-w,cx+w,cy-w,cy+w)); a.plot(cx,cy,'k+',ms=6,mew=1)
    a.set_xticks([]); a.set_yticks([]); a.set_title(f'({lab[j]}) $n={n}$ ($n\\equiv{n%4}$)',loc='left',pad=2)
    for r,pl in ((1,'I'),(2,'II')):
        D=MAPS[pl][n]; V=D['V']; cx,cy=D['centre']
        a=ax[r][j]; Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in Cc)/4; Y=sum(p[1] for p in Cc)/4; W=9
        for c,val in V.items():
            px,py=pcen(c)
            if abs(px-X)<=W+1 and abs(py-Y)<=W+1: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=CMAPV(nv(val)),ec='#666',lw=.25))
        a.plot(X,Y,'k+',ms=6,mew=1); a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
        a.set_title(f'({lab[4*r+j]}) {pattern(n,pl)}',loc='left',pad=2)
ax[0][0].set_ylabel('chessboard, I',fontsize=8); ax[1][0].set_ylabel('Shakti, placement I',fontsize=8); ax[2][0].set_ylabel('Shakti, placement II',fontsize=8)
for r in range(3):
    small_cbar(fig,plt.cm.ScalarMappable(norm=nv,cmap=CMAPV),ax[r][3],label='$\\langle v\\rangle$',size=('6%','80%'),pad=0.04,ticks=[-.5,0,.5])
save(fig,'fig4_exactmaps')
