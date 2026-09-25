import sys, math, random
sys.argv=['x','5','0','0']
src=open('ring_shift.py').read().split("if __name__")[0]
base=open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0]
def build(SX,SY,seed):
    g={}
    code=base.replace("n=5","n=5").replace("random.seed(7)",f"random.seed({seed})")
    # shift the diamond before tiling
    code=code.replace("assert len(cells)==2*n*(n+1)",f"assert len(cells)==2*n*(n+1)\ncells={{(c[0]+{SX},c[1]+{SY}) for c in cells}}")
    exec(code,g)
    exec(src.split("import networkx as nx\n",1)[1].replace("SX,SY=int(sys.argv[2]),int(sys.argv[3])\ncells={(c[0]+SX,c[1]+SY) for c in cells}\n",""),g) if False else None
    exec("import networkx as nx\n"+src.split("cells={(c[0]+SX,c[1]+SY) for c in cells}\n",1)[1],g)
    short={g['short_cell'](d) for d in g['doms']}
    G,bd,ring=g['graph'](short,g['touches_perimeter'])
    return g,G,bd,ring
import networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
def draw(ax,g,G,ring,title):
    rs,inv,outline,doms=g['rs'],g['inv'],g['outline'],g['doms']
    for d in doms:
        sqs=[]
        for c in d:
            r=list(rs[inv[c]]); sqs+=r
            sc=g['short_cell'](d)
            ax.add_patch(Polygon(outline(r),closed=True,fc='#c9c9c9' if c==sc else '#f1f1f1',ec='white',lw=.5))
        ax.add_patch(Polygon(outline(sqs),closed=True,fill=False,ec='#888',lw=1))
    P={c:g['pcen'](c) for c in G}
    bip=nx.is_bipartite(G)
    odd=[]
    if not bip:
        # shortest odd cycle through BFS layering
        best=None
        for s in G:
            dist={s:0};par={s:None};Q=[s]
            for u in Q:
                for v in G[u]:
                    if v not in dist: dist[v]=dist[u]+1;par[v]=u;Q.append(v)
                    elif dist[v]==dist[u] and (best is None or 2*dist[u]+1<len(best)):
                        pu=[u];pv=[v]
                        while par[pu[-1]] is not None: pu.append(par[pu[-1]])
                        while par[pv[-1]] is not None: pv.append(par[pv[-1]])
                        cyc=pu[::-1]+pv  # may share prefix; trim
                        a=pu[::-1];b=pv[::-1];k=0
                        while k<min(len(a),len(b)) and a[k]==b[k]: k+=1
                        cyc=a[k-1:]+b[k:][::-1]
                        if len(cyc)%2==1 and (best is None or len(cyc)<len(best)): best=cyc
        odd=best
    oddE={frozenset(e) for e in zip(odd,odd[1:]+odd[:1])} if odd else set()
    for u,v,dat in G.edges(data=True):
        e=frozenset((u,v)); red=e in oddE
        ax.plot(*zip(P[u],P[v]),color='#d62728' if red else ('#1f77b4' if dat['kind']=='ring' else 'k'),
                ls='--' if dat['kind']=='ring' else '-',lw=3.2 if red else 2,zorder=7 if red else 6)
    side=nx.bipartite.color(G) if bip else {}
    for c in G:
        fc={0:'#f2a93b',1:'#7b4fa0'}.get(side.get(c),'white')
        ax.plot(*P[c],'s' if g['black'](c) else 'o',ms=11,mfc=fc,mec='k',mew=1.5,zorder=8)
    ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=11.5)
    return odd
fig,(A,B)=plt.subplots(1,2,figsize=(15,8))
g,G,bd,ring=build(0,0,7)
same=sum(g['black'](a)==g['black'](b) for a,b in ring)
draw(A,g,G,ring,f'Placement I (the n=5 example): bipartite\n{len(ring)} ring links, {same} of them join parallel plaquettes')
# placement II: find a random tiling that is not bipartite
for seed in range(1,200):
    g2,G2,bd2,ring2=build(1,0,seed)
    if not nx.is_bipartite(G2): break
odd=draw(B,g2,G2,ring2,'')
B.set_title(f'Placement II (diamond shifted by one cell): NOT bipartite\nshortest odd cycle (red), length {len(odd)}',fontsize=11.5)
h=[Line2D([],[],color='k',lw=2,label='interior link (perpendicular contact)'),
   Line2D([],[],color='#1f77b4',ls='--',lw=2,label='added boundary-ring link'),
   Line2D([],[],marker='s',ls='',mfc='w',mec='k',ms=9,label='short arm on horizontal plaquette'),
   Line2D([],[],marker='o',ls='',mfc='w',mec='k',ms=9,label='short arm on vertical plaquette'),
   Line2D([],[],marker='o',ls='',mfc='#f2a93b',mec='k',ms=9,label='bipartition side A'),
   Line2D([],[],marker='o',ls='',mfc='#7b4fa0',mec='k',ms=9,label='bipartition side B')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=10)
plt.tight_layout(rect=(0,0.08,1,0.97))
for ext in ('png','pdf'): plt.savefig(f'outputs/aztec_n5_ring_bipartite.{ext}',dpi=200)
print('seed',seed)
