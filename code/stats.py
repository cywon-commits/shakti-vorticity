exec(open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0].replace("n=5","n=int(__import__('os').environ.get('NN','5'))").replace("random.seed(7)","random.seed(int(__import__('os').environ.get('SEED','7')))").replace("range(20000)","range(int(__import__('os').environ.get('IT','20000')))"))
import networkx as nx
short={short_cell(d) for d in doms}
G=nx.Graph(); G.add_nodes_from(short)
for c in short:
    for dx,dy in [(1,0),(0,1)]:
        if (c[0]+dx,c[1]+dy) in short: G.add_edge(c,(c[0]+dx,c[1]+dy))
print(max(d for _,d in G.degree()), len(nx.cycle_basis(G)), nx.number_connected_components(G))
