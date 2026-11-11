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

def simplify(poly):
    out=[]
    m=len(poly)
    for k in range(m):
        a,b,c=poly[k-1],poly[k],poly[(k+1)%m]
        if (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0])!=0: out.append(b)
    return out
def inset(poly,d):
    """move every edge of a CCW rectilinear loop inward by d"""
    p=simplify(poly); m=len(p); q=[]
    for k in range(m):
        a,b,c=p[k-1],p[k],p[(k+1)%m]
        e1=(b[0]-a[0],b[1]-a[1]); e2=(c[0]-b[0],c[1]-b[1])
        u=lambda e:(e[0]/max(abs(e[0]),abs(e[1])),e[1]/max(abs(e[0]),abs(e[1])))
        n1=(-u(e1)[1],u(e1)[0]); n2=(-u(e2)[1],u(e2)[0])
        q.append((b[0]+d*(n1[0]+n2[0]),b[1]+d*(n1[1]+n2[1])))
    return q
GAP=0.13
fig=plt.figure(figsize=(TWO,4.45))
gs=fig.add_gridspec(1,3,width_ratios=[1,1.02,1.02],wspace=0.06,left=0.005,right=0.995,top=0.955,bottom=0.405)
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
        dd=[(c[0]+sh[0],c[1]+sh[1]) for c in d]; t=dtype(d); sq=[]   # hue = type in (a): same hue = same domino in all panels
        rs=[]
        for c0,c in zip(d,dd):
            r=list(rect(c)); sq+=r; used.append(r); rs.append((r,c))
        shape=Polygon(inset(outline(sq),GAP),closed=True,fill=False,ec='k',lw=1.0,zorder=6)
        clip=Polygon(inset(outline(sq),GAP),closed=True,transform=B.transData)
        for r,c in rs:
            P=B.add_patch(Polygon(outline(r),closed=True,fc=shade(typecol[t],black(c)),ec='none',zorder=3))
            P.set_clip_path(clip)
        # thin white line where the two plaquettes of the tetromino meet
        r0=set(map(tuple,rs[0][0])); r1=set(map(tuple,rs[1][0]))
        for (x,y) in r0:
            for (dx,dy) in ((1,0),(-1,0),(0,1),(0,-1)):
                if (x+dx,y+dy) in r1:
                    if dx: xx=x+max(dx,0); L=B.plot([xx,xx],[y+GAP*0,y+1],color='white',lw=.8,zorder=4)
                    else: yy=y+max(dy,0); L=B.plot([x,x+1],[yy,yy],color='white',lw=.8,zorder=4)
                    L[0].set_clip_path(clip)
        B.add_patch(shape)
        el=next(p for p in sq if (lambda nb: len(nb)==2 and nb[0][0]*nb[1][0]+nb[0][1]*nb[1][1]==0)([(q[0]-p[0],q[1]-p[1]) for q in sq if abs(p[0]-q[0])+abs(p[1]-q[1])==1]))
        B.text(el[0]+.5,el[1]+.5,str(k),ha='center',va='center',fontsize=5.5,zorder=5)
    edges=set()
    for r in used:
        for (x,y) in r:
            for e in [((x,y),(x+1,y)),((x,y+1),(x+1,y+1)),((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]:
                if is_island(e): edges.add(e)
    for (a,b) in edges:
        if a[0]==b[0]: B.plot([a[0],a[0]],[a[1]+.18,b[1]-.18],color='#6b6b6b',lw=1.9,solid_capstyle='butt',zorder=2)
        else: B.plot([a[0]+.18,b[0]-.18],[a[1],a[1]],color='#6b6b6b',lw=1.9,solid_capstyle='butt',zorder=2)
    xs=[p[0] for r in used for p in r]; ys=[p[1] for r in used for p in r]
    B.set_xlim(min(xs)-.3,max(xs)+1.3); B.set_ylim(min(ys)-.3,max(ys)+1.3); B.set_aspect('equal'); B.axis('off')
    B.set_title(f'Shakti lattice, placement {pl}',pad=2); label_in(B,'(b)' if idx==0 else '(c)',x=0.0,y=1.0)
h=[Patch(fc=shade(typecol[t],True),ec='none',label=lab) for t,lab in (('E','E'),('W','W'),('N','N'),('S','S'))]
fig.legend(handles=h+[Patch(fc='#999',label='dark: black cell = horizontal plaquette'),Patch(fc='#ddd',label='light: white cell = vertical plaquette'),
                     plt.Line2D([],[],color='#6b6b6b',lw=2,label='island')],loc='lower center',ncol=7,fontsize=6.5,handlelength=1.2,columnspacing=0.9,bbox_to_anchor=(0.5,0.335))

# ---- (d) geometry of the map: chessboard square -> Shakti diamond, chessboard diamond -> Shakti square
def chess_ax(a,cells):
    for c in cells: a.add_patch(Polygon(outline([c]),closed=True,fc='#7d8fb3' if black(c) else '#c9d3e6',ec='white',lw=.3))
    a.add_patch(Polygon(outline(list(cells)),closed=True,fill=False,ec='k',lw=.8))
    xs=[c[0] for c in cells]; ys=[c[1] for c in cells]
    cx=(min(xs)+max(xs)+1)/2; cy=(min(ys)+max(ys)+1)/2; a.set_xlim(cx-HR,cx+HR); a.set_ylim(cy-HR,cy+HR); a.set_aspect('equal'); a.axis('off')
def shakti_ax(a,cells):
    sq=[]
    for c in cells:
        r=list(rect(c)); sq+=r
        a.add_patch(Polygon(outline(r),closed=True,fc='#7d8fb3' if black(c) else '#c9d3e6',ec='white',lw=.3))
    a.add_patch(Polygon(outline(sq),closed=True,fill=False,ec='k',lw=.8))
    xs=[p[0] for p in sq]; ys=[p[1] for p in sq]
    cx=(min(xs)+max(xs)+1)/2; cy=(min(ys)+max(ys)+1)/2; a.set_xlim(cx-HR,cx+HR); a.set_ylim(cy-HR,cy+HR); a.set_aspect('equal'); a.axis('off')
HR=8.6   # same data range in all four panels: lengths are to scale
SQ={(i,j) for i in range(-4,4) for j in range(-4,4)}
AZ={(i+1,j) for i in range(-4,4) for j in range(-4,4) if abs(i+.5)+abs(j+.5)<=4}
gd=fig.add_gridspec(1,7,width_ratios=[1,.42,1,.45,1,.42,1],wspace=0.0,left=0.06,right=0.94,top=0.25,bottom=0.01)
P=[fig.add_subplot(gd[k]) for k in (0,2,4,6)]
chess_ax(P[0],SQ); shakti_ax(P[1],SQ); chess_ax(P[2],AZ); shakti_ax(P[3],AZ)
fig.canvas.draw()
for a,t in zip(P,('chessboard square, $8\\times8$','Shakti: rotated by $45^\\circ$','chessboard diamond (Aztec, $n=4$)','Shakti: axis-aligned')):
    bb=a.get_position(); fig.text((bb.x0+bb.x1)/2,0.262,t,ha='center',va='bottom',fontsize=7)
for k in (1,5):
    bb=fig.add_subplot(gd[k]); bb.axis('off')
    bb.annotate('',xy=(0.95,0.5),xytext=(0.05,0.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='-|>',lw=1.1,color='k'))
    bb.text(0.5,0.62,'$45^\\circ$, $\\times\\sqrt{2}$',ha='center',va='bottom',fontsize=6.5,transform=bb.transAxes)
fig.text(0.005,0.262,'(d)',ha='left',va='bottom',fontsize=9,fontweight='bold')
save(fig,'fig1_mapping')
