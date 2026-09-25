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
M=pickle.load(open('exact/maps.pkl','rb'))
fig,ax=plt.subplots(2,4,figsize=(TWO,3.55),gridspec_kw=dict(wspace=0.05,hspace=0.12,left=0.02,right=0.93,top=0.95,bottom=0.0))
for j,n in enumerate([36,37,38,39]):
    D=M[n]; V=D['V']; cx,cy=D['centre']; w=8
    img=np.full((2*w,2*w),np.nan)
    for i in range(-w,w):
        for k in range(-w,w):
            if (cx+i,cy+k) in V: img[k+w,i+w]=V[(cx+i,cy+k)]
    a=ax[0][j]; a.imshow(img,origin='lower',cmap=CMAPV,norm=nv,extent=(cx-w,cx+w,cy-w,cy+w)); a.plot(cx,cy,'k+',ms=6,mew=1)
    a.set_xticks([]); a.set_yticks([]); a.set_title(f'({"abcd"[j]}) $n={n}$ ($n\\equiv{n%4}$)',loc='left',pad=2)
    a=ax[1][j]; Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in Cc)/4; Y=sum(p[1] for p in Cc)/4; W=9
    for c,val in V.items():
        px,py=pcen(c)
        if abs(px-X)<=W+1 and abs(py-Y)<=W+1: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=CMAPV(nv(val)),ec='#666',lw=.25))
    a.plot(X,Y,'k+',ms=6,mew=1); a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'({"efgh"[j]}) Shakti view',loc='left',pad=2)
ax[0][0].set_ylabel('chessboard',fontsize=8); ax[1][0].set_ylabel('Shakti lattice',fontsize=8)
small_cbar(fig,plt.cm.ScalarMappable(norm=nv,cmap=CMAPV),ax[0][3],label='$\\langle v\\rangle$',size=('6%','80%'),pad=0.04,ticks=[-.5,0,.5])
small_cbar(fig,plt.cm.ScalarMappable(norm=nv,cmap=CMAPV),ax[1][3],label='$\\langle v\\rangle$',size=('6%','80%'),pad=0.04,ticks=[-.5,0,.5])
save(fig,'fig4_exactmaps')
