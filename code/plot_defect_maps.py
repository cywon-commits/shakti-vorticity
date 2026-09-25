import pickle, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect, pcen
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
cmapr=LinearSegmentedColormap.from_list('rel',['#6a1b9a','#b39ddb','#f4f4f4','#8fd19e','#1b7a3a']); nr=Normalize(-1.2,1.2)
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
def frame(D):
    cx,cy=D['centre']; Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]
    return sum(p[0] for p in Cc)/4, sum(p[1] for p in Cc)/4
def panel(a,D,vals,cmap,norm,mask=None,defect=None,title=''):
    X,Y=frame(D); Wd=9.5
    for c,val in vals.items():
        px,py=pcen(c)
        if abs(px-X)>Wd or abs(py-Y)>Wd: continue
        fc='#cfcfcf' if (mask is not None and c in mask) else cmap(norm(val))
        a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=fc,ec='#666',lw=.35))
    if defect=='O':
        for O in D['Os']:
            for c in O: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#555',ec='k',lw=.8))
            p=[pcen(c) for c in O]; a.text((p[0][0]+p[1][0])/2,(p[0][1]+p[1][1])/2,'O',color='w',ha='center',va='center',fontsize=9,fontweight='bold')
    if defect=='mono':
        A,B=D['AB']
        for c in (A,B):
            a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#555',ec='k',lw=.8)); a.plot(*pcen(c),'*',color='w',ms=10)
        pts=[pcen(A)]+[((pcen((x-1,A[1]))[0]+pcen((x,A[1]))[0])/2,(pcen((x-1,A[1]))[1]+pcen((x,A[1]))[1])/2) for x in range(A[0]+1,B[0]+1)]+[pcen(B)]
        a.plot([p[0] for p in pts],[p[1] for p in pts],'k--',lw=1.2)
    a.set_xlim(X-Wd,X+Wd); a.set_ylim(Y-Wd,Y+Wd); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title(title,fontsize=10.5)
for pl in ('I','II'):
    M=pickle.load(open(f'exact/defect_maps_{pl}.pkl','rb'))
    fig,ax=plt.subplots(5,4,figsize=(17,22))
    for j,n in enumerate(sorted(M)):
        D=M[n]; V0=D['none']
        panel(ax[0][j],D,V0,cmapv,nv,title=f'n = {n} (n ≡ {n%4}): no defect')
        panel(ax[1][j],D,D['O'],cmapv,nv,defect='O',title='O pair (full vortices)')
        panel(ax[2][j],D,D['mono'],cmapv,nv,defect='mono',title='monomer pair (half vortices)')
        for row,key,df in ((3,'O','O'),(4,'mono','mono')):
            rel={}; mask=set()
            for c,val in D[key].items():
                v0=V0.get(c,0)
                if abs(v0)<0.04: mask.add(c); rel[c]=0
                else: rel[c]=val/v0
            panel(ax[row][j],D,rel,cmapr,nr,mask=mask,defect=df,title=('O pair' if key=='O' else 'monomer pair')+': ⟨v⟩ / ⟨v⟩₀')
    fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax[:3].ravel().tolist(),shrink=.5,pad=.02).set_label('⟨v⟩  (orange +1 = counter-clockwise, blue −1 = clockwise)')
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=nr,cmap=cmapr),ax=ax[3:].ravel().tolist(),shrink=.8,pad=.02); cb.set_label('⟨v⟩ / ⟨v⟩₀:  +1 unchanged, 0 erased, −1 reversed  (grey: |⟨v⟩₀| < 0.04)')
    fig.suptitle(f'Exact ⟨v⟩ with fixed defects, Aztec placement {pl} (Kasteleyn; defects 7 cells apart around the centre)',fontsize=14,y=.995)
    plt.savefig(f'outputs/defect_maps_{pl}.png',dpi=110); plt.savefig(f'outputs/defect_maps_{pl}.pdf')
    print('saved',pl)
