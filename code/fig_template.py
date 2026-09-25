import io, contextlib, math
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrow
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
def arrow(ax,a,b,color,lw=1.2,double=False,fill=True):
    (x0,y0),(x1,y1)=a,b; dx,dy=x1-x0,y1-y0; m=0.18
    sx,sy=x0+m*dx,y0+m*dy; ex,ey=x1-m*dx,y1-m*dy
    if double:
        ax.annotate('',xy=(ex,ey),xytext=(sx,sy),arrowprops=dict(arrowstyle='<|-|>',color=color,lw=2.2,mutation_scale=14))
    else:
        ax.annotate('',xy=(ex,ey),xytext=(sx,sy),arrowprops=dict(arrowstyle='-|>',color=color,lw=2.6 if fill else 1.6,mutation_scale=16,
                    fc=color if fill else 'white',ec='k' if not fill else color))
def draw(ax,mirror,v,title):
    f=(lambda p:(2-p[0],p[1])) if mirror else (lambda p:p)
    s=-1 if mirror else 1
    sq=lambda c:(1-c[0],c[1]) if mirror else c
    short=[sq(c) for c in T_SHORT]; long_=[sq(c) for c in T_SQ if c not in T_SHORT]
    col_s,col_l=(('#e8743a','#f6c6a8') if v>0 else ('#2f6db5','#aac6ea'))
    ax.add_patch(Polygon(outline(short),closed=True,fc=col_s,ec='k',lw=1.5))
    ax.add_patch(Polygon(outline(long_),closed=True,fc=col_l,ec='k',lw=1.5))
    k=s*v
    for i,(a,b) in enumerate(T_PER):
        a2,b2=(f(a),f(b)) if k>0 else (f(b),f(a))
        arrow(ax,a2,b2,'#c00000' if i==8 else '#2e9d4f')
    for (a,b) in T_INT+T_EXT:
        a2,b2=(f(a),f(b)) if k>0 else (f(b),f(a)); arrow(ax,a2,b2,'#666',fill=False)
    for (a,b) in T_FREE: arrow(ax,f(a),f(b),'k',double=True)
    ax.set_xlim(-1.4,3.4); ax.set_ylim(-1.3,4.3); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=12)
fig,ax=plt.subplots(1,4,figsize=(18,5.5))
draw(ax[0],False,+1,'L, v = +1 (counter-clockwise)'); draw(ax[1],False,-1,'L, v = −1 (clockwise)')
draw(ax[2],True,+1,'J, v = +1'); draw(ax[3],True,-1,'J, v = −1')
from matplotlib.lines import Line2D
h=[Line2D([],[],color='#2e9d4f',lw=3,label='perimeter island along the circulation'),Line2D([],[],color='#c00000',lw=3,label='defect island (against the circulation)'),
   Line2D([],[],color='#666',lw=1.6,label='internal / neighbouring island fixed by the template'),Line2D([],[],color='k',lw=2.2,label='free island'),
   Polygon([[0,0]],fc='#e8743a',label='strong (short) arm'),Polygon([[0,0]],fc='#f6c6a8',label='weak (long) arm')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=10)
plt.tight_layout(rect=(0,0.12,1,0.93)); plt.savefig('outputs/fig_template.png',dpi=180); plt.savefig('outputs/fig_template.pdf')
# ---- flip lemma figure: 2x2 chess blocks of both types, both fillings, in the Shakti view
fig,ax=plt.subplots(2,2,figsize=(10,10))
for row,(x,y,typ) in enumerate([(0,0,'x+y even: opposite v in both fillings (ΣV = 0)'),(1,0,'x+y odd: equal v in both fillings (ΣV = ±2, same sign)')]):
    for col,fill in enumerate(('h','v')):
        a=ax[row][col]
        if fill=='h': ds=[((x,y),(x+1,y)),((x,y+1),(x+1,y+1))]
        else: ds=[((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]
        vs=(1,-1) if row==0 else (1,1)
        allp=[]
        for kk,(d,v) in enumerate(zip(ds,vs)):
            sq=[p for c in d for p in rect(c)]; allp+=sq
            fc=('#e8743a' if v>0 else '#2f6db5') if kk==0 else ('#f0a070' if v>0 else '#6f9bd6')
            a.add_patch(Polygon(outline(sq),closed=True,fc=fc,ec='k',lw=2.5,hatch='' if kk==0 else '//'))
            for c in d: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fill=False,ec='white',lw=.8))
            cx=sum(p[0]+.5 for p in sq)/4; cy=sum(p[1]+.5 for p in sq)/4
            a.text(cx,cy,'↺' if v>0 else '↻',ha='center',va='center',fontsize=18,color='white')
        xs=[p[0] for p in allp]; ys=[p[1] for p in allp]
        a.set_xlim(min(xs)-.5,max(xs)+1.5); a.set_ylim(min(ys)-.5,max(ys)+1.5); a.set_aspect('equal'); a.axis('off')
        a.set_title(f'{"horizontal" if fill=="h" else "vertical"} dominoes', fontsize=11)
    ax[row][0].text(-0.15,0.5,typ,transform=ax[row][0].transAxes,rotation=90,va='center',ha='center',fontsize=11)
fig.suptitle('Lemma 3: a 2×2 flip in the Shakti view (colours: one representative sign choice)',fontsize=12)
plt.tight_layout(rect=(0.03,0,1,0.96)); plt.savefig('outputs/fig_flip_lemma.png',dpi=170); plt.savefig('outputs/fig_flip_lemma.pdf')
