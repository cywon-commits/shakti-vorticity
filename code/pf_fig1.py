import random, io, contextlib
from pubstyle import *
from matplotlib.patches import Polygon, Patch
from matplotlib.colors import to_rgb
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
from vort import is_island
n=5
base={(i,j) for i in range(-n,n) for j in range(-n,n) if abs(i+.5)+abs(j+.5)<=n}
dom={}
for (i,j) in sorted(base):
    if (i,j) in dom: continue
    dom[(i,j)]=(i+1,j); dom[(i+1,j)]=(i,j)
random.seed(7); cl=sorted(base)
for _ in range(20000):
    i,j=random.choice(cl); sq=[(i,j),(i+1,j),(i,j+1),(i+1,j+1)]
    if not all(c in base for c in sq): continue
    if dom[(i,j)]==(i+1,j) and dom[(i,j+1)]==(i+1,j+1):
        dom[(i,j)]=(i,j+1);dom[(i,j+1)]=(i,j);dom[(i+1,j)]=(i+1,j+1);dom[(i+1,j+1)]=(i+1,j)
    elif dom[(i,j)]==(i,j+1) and dom[(i+1,j)]==(i+1,j+1):
        dom[(i,j)]=(i+1,j);dom[(i+1,j)]=(i,j);dom[(i,j+1)]=(i+1,j+1);dom[(i+1,j+1)]=(i,j+1)
doms0=sorted({tuple(sorted((c,dom[c]))) for c in base})
black=lambda c:(c[0]+c[1])%2==0
typecol={'E':CAT[1],'W':CAT[0],'N':CAT[3],'S':CAT[2]}
def dtype(d):
    b=d[0] if black(d[0]) else d[1]; w=d[1] if b==d[0] else d[0]
    return {(1,0):'E',(-1,0):'W',(0,1):'N',(0,-1):'S'}[(w[0]-b[0],w[1]-b[1])]
def shade(c,bl):
    r=to_rgb(c); return tuple(0.8*x for x in r) if bl else tuple(0.45*x+0.55 for x in r)
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
fig=plt.figure(figsize=(TWO,2.75))
gs=fig.add_gridspec(1,3,width_ratios=[1,1.02,1.02],wspace=0.06,left=0.005,right=0.995,top=0.93,bottom=0.12)
A=fig.add_subplot(gs[0])
for k,d in enumerate(doms0,1):
    t=dtype(d)
    for c in d: A.add_patch(Polygon(outline([c]),closed=True,fc=shade(typecol[t],black(c)),ec='white',lw=.4))
    A.add_patch(Polygon(outline(list(d)),closed=True,fill=False,ec='k',lw=.9))
    A.text(sum(c[0] for c in d)/2+.5,sum(c[1] for c in d)/2+.5,str(k),ha='center',va='center',fontsize=5.5)
A.set_xlim(-n-.3,n+.3); A.set_ylim(-n-.3,n+.3); A.set_aspect('equal'); A.axis('off')
A.set_title('Aztec diamond, $n=5$',pad=2); label_in(A,'(a)',x=0.0,y=1.0)
for idx,(pl,sh) in enumerate((('I',((n+1)%2,0)),('II',(n%2,0)))):
    B=fig.add_subplot(gs[idx+1]); used=[]
    for k,d in enumerate(doms0,1):
        dd=[(c[0]+sh[0],c[1]+sh[1]) for c in d]; t=dtype(tuple(dd)); sq=[]
        for c0,c in zip(d,dd):
            r=list(rect(c)); sq+=r; used.append(r)
            B.add_patch(Polygon(outline(r),closed=True,fc=shade(typecol[t],black(c)),ec='none'))
        B.add_patch(Polygon(outline(sq),closed=True,fill=False,ec='k',lw=.9,zorder=3))
        el=next(p for p in sq if (lambda nb: len(nb)==2 and nb[0][0]*nb[1][0]+nb[0][1]*nb[1][1]==0)([(q[0]-p[0],q[1]-p[1]) for q in sq if abs(p[0]-q[0])+abs(p[1]-q[1])==1]))
        B.text(el[0]+.5,el[1]+.5,str(k),ha='center',va='center',fontsize=5.5,zorder=5)
    edges=set()
    for r in used:
        for (x,y) in r:
            for e in [((x,y),(x+1,y)),((x,y+1),(x+1,y+1)),((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]:
                if is_island(e): edges.add(e)
    for (a,b) in edges:
        if a[0]==b[0]: B.plot([a[0],a[0]],[a[1]+.18,b[1]-.18],color='#555',lw=2.0,solid_capstyle='round',zorder=4)
        else: B.plot([a[0]+.18,b[0]-.18],[a[1],a[1]],color='#555',lw=2.0,solid_capstyle='round',zorder=4)
    xs=[p[0] for r in used for p in r]; ys=[p[1] for r in used for p in r]
    B.set_xlim(min(xs)-.3,max(xs)+1.3); B.set_ylim(min(ys)-.3,max(ys)+1.3); B.set_aspect('equal'); B.axis('off')
    B.set_title(f'Shakti lattice, placement {pl}',pad=2); label_in(B,'(b)' if idx==0 else '(c)',x=0.0,y=1.0)
h=[Patch(fc=shade(typecol[t],True),ec='none',label=lab) for t,lab in (('E','E'),('W','W'),('N','N'),('S','S'))]
fig.legend(handles=h+[Patch(fc='#999',label='dark: black cell = horizontal plaquette'),Patch(fc='#ddd',label='light: white cell = vertical plaquette'),
                     plt.Line2D([],[],color='#555',lw=2,label='island')],loc='lower center',ncol=7,fontsize=6.5,handlelength=1.2,columnspacing=0.9,bbox_to_anchor=(0.5,0.0))
save(fig,'fig1_mapping')
