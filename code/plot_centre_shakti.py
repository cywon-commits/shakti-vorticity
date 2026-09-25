import pickle, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect, pcen, black
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
W=9
fig,ax=plt.subplots(2,2,figsize=(13,13)); ax=ax.ravel()
for a,n in zip(ax,[36,37,38,39]):
    D=pickle.load(open(f'fast/f_{n}_1500.pkl','rb')); x0,y0=D['x0'],D['y0']; cx,cy=D['center']; V=D['vmean']
    # Shakti position of the chessboard centre vertex = mean of the four central cells' plaquette centres
    C=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]
    X=sum(p[0] for p in C)/4; Y=sum(p[1] for p in C)/4
    for c in D['cells']:
        px,py=pcen(c)
        if abs(px-X)<=W and abs(py-Y)<=W:
            val=V[c[0]-x0,c[1]-y0]
            a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=cmapv(nv(val)),ec='#555',lw=.5))
    a.plot(X,Y,'k+',ms=14,mew=2)
    a.set_xlim(X-W,X+W); a.set_ylim(Y-W,Y+W); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'n = {n} (n≡{n%4} mod 4):  ⟨v⟩ per plaquette',fontsize=12)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax.tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩   (blue −1 clockwise, orange +1 counter-clockwise)')
fig.suptitle('Centre of the Aztec diamond on the Shakti lattice\n(+ = centre, 1500 uniform tilings per n; same data as the chessboard panels)',fontsize=13)
plt.savefig('outputs/vorticity_centre_shakti.png',dpi=170); plt.savefig('outputs/vorticity_centre_shakti.pdf')
