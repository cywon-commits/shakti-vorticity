exec(open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0])
import networkx as nx
short={short_cell(d):k for k,d in enumerate(doms,1)}          # chess cell -> domino number
# edge: two short-arm plaquettes that are perpendicular edge-neighbours (= orthogonal chess neighbours).
# parallel neighbours (same block, same orientation) are NOT connected.
G=nx.Graph(); G.add_nodes_from(short.values())
for c,k in short.items():
    for dx,dy in [(1,0),(0,1)]:
        c2=(c[0]+dx,c[1]+dy)
        if c2 in short:
            assert isL(rs[inv[c]],rs[inv[c2]]) and black(c)!=black(c2)
            G.add_edge(k,short[c2])
# sanity: parallel short-arm pairs that were excluded
par=[(short[c],short[c2]) for c in short for c2 in short if c<c2 and inv[c][:2]==inv[c2][:2]]
comps=sorted(nx.connected_components(G),key=len,reverse=True)
print('edges',G.number_of_edges(),'components',[sorted(x) for x in comps]); print('excluded parallel pairs',par)
print('degrees',dict(sorted(G.degree())))

def pcen(c):
    r=rs[inv[c]]; return (sum(p[0] for p in r)/2+.5, sum(p[1] for p in r)/2+.5)
fig,(A,B)=plt.subplots(1,2,figsize=(15,7.8))
for ax,shak in [(A,False),(B,True)]:
    for k,d in enumerate(doms,1):
        t=dtype(d); sc=short_cell(d); sqs=[]
        for c in d:
            poly=outline(list(rs[inv[c]])) if shak else outline([c]); sqs+= list(rs[inv[c]]) if shak else [c]
            ax.add_patch(Polygon(poly,closed=True,fc=col[t] if c==sc else pale(col[t]),ec='white',lw=.5,alpha=.55 if c==sc else 1))
        ax.add_patch(Polygon(outline(sqs),closed=True,fill=False,ec='#333',lw=1.3))
    pos={k:(pcen(c) if shak else (c[0]+.5,c[1]+.5)) for c,k in short.items()}
    for u,v in par:
        ax.plot(*zip(pos[u],pos[v]),ls=':',color='#999',lw=1.5,zorder=6)
    for u,v in G.edges(): ax.plot(*zip(pos[u],pos[v]),color='k',lw=2.4,zorder=7)
    for c,k in short.items():
        ax.plot(*pos[k],'o',ms=15,mfc='k' if black(c) else 'white',mec='k',mew=1.8,zorder=8)
        ax.text(*pos[k],str(k),ha='center',va='center',fontsize=7.5,color='white' if black(c) else 'k',fontweight='bold',zorder=9)
    ax.set_aspect('equal');ax.axis('off')
A.set_xlim(-n-.5,n+.5);A.set_ylim(-n-.5,n+.5)
xs=[p[0] for c in cells for p in rs[inv[c]]]; ys=[p[1] for c in cells for p in rs[inv[c]]]
B.set_xlim(min(xs)-.5,max(xs)+1.5);B.set_ylim(min(ys)-.5,max(ys)+1.5)
A.set_title('Chessboard view: short-arm cells and their links',fontsize=12)
B.set_title('Shakti view: graph of short-arm plaquettes\n(link = perpendicular plaquettes sharing an edge)',fontsize=12)
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',mfc='k',mec='k',ms=10,label='short arm on horizontal plaquette (black cell)'),
   Line2D([],[],marker='o',ls='',mfc='w',mec='k',ms=10,label='short arm on vertical plaquette (white cell)'),
   Line2D([],[],color='k',lw=2.4,label=f'link: perpendicular contact ({G.number_of_edges()} links)'),
   Line2D([],[],color='#999',ls=':',lw=1.5,label=f'parallel contact, not linked ({len(par)})')]
fig.legend(handles=h,loc='lower center',ncol=2,frameon=False,fontsize=10)
plt.tight_layout(rect=(0,0.08,1,0.97))
for ext in ('png','pdf'): plt.savefig(f'outputs/aztec_n5_shortarm_graph.{ext}',dpi=200)
