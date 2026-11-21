import io, contextlib
from pubstyle import *
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
with contextlib.redirect_stdout(io.StringIO()):
    from vort import T_PER,T_INT,T_EXT,T_FREE,T_SQ,T_SHORT
    from core import rect
def outline(sqs):
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; s=next(iter(nxt)); poly=[s]; p=nxt[s]
    while p!=s: poly.append(p); p=nxt[p]
    return poly
import numpy as np
import matplotlib.patheffects as pe
CAPF={'per':'#dfeee4','def':'#f6d9d9','fix':'#ececec','free':'#ffffff'}
CAPA={'per':'#2e9d4f','def':'#c00000','fix':'#6b6b6b','free':'#333333'}
def island(ax,a,b,role):
    a=np.array(a,float); b=np.array(b,float); d=b-a; m=0.17
    s,e=a+m*d,b-m*d; u=d/np.hypot(*d)
    ax.plot(*zip(s,e),color='#555',lw=6.2,solid_capstyle='round',zorder=4)
    ax.plot(*zip(s,e),color=CAPF[role],lw=5.1,solid_capstyle='round',zorder=4)
    s2,e2=s-0.02*u,e+0.02*u
    if role=='free':
        ax.annotate('',xy=e2,xytext=s2,arrowprops=dict(arrowstyle='<|-|>,head_width=0.18,head_length=0.32',color=CAPA[role],lw=0.8,shrinkA=0,shrinkB=0),zorder=5)
    else:
        ax.annotate('',xy=e2,xytext=s2,arrowprops=dict(arrowstyle='-|>,head_width=0.2,head_length=0.38',color=CAPA[role],lw=1.3,shrinkA=0,shrinkB=0),zorder=5)
def draw(ax,mirror,v,title):
    f=(lambda p:(2-p[0],p[1])) if mirror else (lambda p:p); s=-1 if mirror else 1
    sq=lambda c:(1-c[0],c[1]) if mirror else c
    short=[sq(c) for c in T_SHORT]; long_=[sq(c) for c in T_SQ if c not in T_SHORT]
    cs,cl=((OR,ORl) if v>0 else (BL,BLl))
    ax.add_patch(Polygon(outline(short),closed=True,fc=cs,ec='k',lw=.6,zorder=1)); ax.add_patch(Polygon(outline(long_),closed=True,fc=cl,ec='k',lw=.6,zorder=1))
    k=s*v
    for i,(a,b) in enumerate(T_PER):
        a2,b2=(f(a),f(b)) if k>0 else (f(b),f(a)); island(ax,a2,b2,'def' if i==8 else 'per')
    for (a,b) in T_INT+T_EXT:
        a2,b2=(f(a),f(b)) if k>0 else (f(b),f(a)); island(ax,a2,b2,'fix')
    for (a,b) in T_FREE: island(ax,f(a),f(b),'free')
    ax.set_xlim(-1.25,3.25); ax.set_ylim(-1.15,4.15); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,loc='left',pad=1)
fig,ax=plt.subplots(1,4,figsize=(TWO,2.35),gridspec_kw=dict(wspace=0.0,left=0,right=1,top=0.93,bottom=0.16))
for a,(m,v,t),lab in zip(ax,[(False,1,'L, $v=+1$'),(False,-1,'L, $v=-1$'),(True,1,'J, $v=+1$'),(True,-1,'J, $v=-1$')],'abcd'):
    draw(a,m,v,f'({lab}) '+t)
cap=lambda r,l: Line2D([],[],color=CAPF[r],lw=5.1,solid_capstyle='round',label=l,path_effects=[pe.Stroke(linewidth=6.2,foreground='#555'),pe.Normal()])
h=[cap('per','perimeter island, moment along the circulation'),cap('def','defect island, moment against it'),
   cap('fix','other island fixed by the template'),cap('free','free island (double arrow)'),
   Patch(fc=OR,label='strong arm'),Patch(fc=ORl,label='weak arm')]
fig.legend(handles=h,loc='lower center',ncol=3,fontsize=7,handlelength=1.6,columnspacing=1.2,bbox_to_anchor=(0.5,-0.01))
save(fig,'fig2_template')
# flip lemma (appendix), single column
fig,ax=plt.subplots(2,2,figsize=(COL,3.3),gridspec_kw=dict(wspace=0.02,hspace=0.12,left=0.1,right=1,top=0.95,bottom=0))
for row,(x,y,typ) in enumerate([(0,0,'$x+y$ even: $\\sum v=0$'),(1,0,'$x+y$ odd: $\\sum v=\\pm2$')]):
    for col,fill in enumerate(('h','v')):
        a=ax[row][col]
        ds=[((x,y),(x+1,y)),((x,y+1),(x+1,y+1))] if fill=='h' else [((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]
        vs=(1,-1) if row==0 else (1,1); allp=[]
        for kk,(d,v) in enumerate(zip(ds,vs)):
            sqs=[p for c in d for p in rect(c)]; allp+=sqs
            fc=(OR if v>0 else BL) if kk==0 else ('#f0a070' if v>0 else '#6f9bd6')
            a.add_patch(Polygon(outline(sqs),closed=True,fc=fc,ec='k',lw=1.2,hatch='' if kk==0 else '////'))
            for c in d: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fill=False,ec='white',lw=.5))
            cx=sum(p[0]+.5 for p in sqs)/4; cy=sum(p[1]+.5 for p in sqs)/4
            a.text(cx,cy,'↺' if v>0 else '↻',ha='center',va='center',fontsize=11,color='white',fontweight='bold')
        xs=[p[0] for p in allp]; ys=[p[1] for p in allp]
        a.set_xlim(min(xs)-.3,max(xs)+1.3); a.set_ylim(min(ys)-.3,max(ys)+1.3); a.set_aspect('equal'); a.axis('off')
        if row==0: a.set_title('horizontal pair' if fill=='h' else 'vertical pair',pad=2)
    ax[row][0].text(-0.08,0.5,typ,transform=ax[row][0].transAxes,rotation=90,va='center',ha='center',fontsize=8)
save(fig,'fig_flip')
