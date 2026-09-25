from mono_ilp import *
from plaq_count import *
M=Model(8); T=M.T
for seed in range(3):
    d=M.tiling(rng=random.Random(seed)); E,sp,labs,el=M.solve(d)
    L,res=plaq_account(T,sp); print(E,Counter(r[0] for r in res))
    # for each domino: the elbow vertex vs plaquettes touching it
    c=d[0]; e=M.elb[(c[0],c[1])]
    touch=[i for i,r in enumerate(T.plaqs) if e in corners(T,r)]
    print(' elbow touches plaquettes',touch,'domino plaquettes',c[:2])
    # angle-based: is the elbow at a corner (90deg) of plaquette or on straight side?
    for i in touch:
        r=T.plaqs[i]; sq=[T.red(q) for q in r]
        n=sum(1 for (x,y) in sq for dx in (0,1) for dy in (0,1) if T.red((x+dx,y+dy))==e)
        print('   plaquette',i,'elbow shared by',n,'unit squares of it (1=corner of rect, 2=middle of long side)')
