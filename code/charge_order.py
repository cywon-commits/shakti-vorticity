from mono_ilp import *
M=Model(8); T=M.T; pats=set()
# enumerate GS via random objective perturbation
import random
for seed in range(12):
    rng=random.Random(seed); cost=dict(COST0)
    prob_doms=M.tiling(rng=rng)
    E,sp,labs,el=M.solve(prob_doms)
    if E!=16: continue
    q=tuple(T.vstate(p,sp)[1] for p in sorted(T.verts) if len(T.inc[p])==3)
    if q[0]<0: q=tuple(-x for x in q)
    pats.add(q)
print('distinct z3 charge patterns (up to global sign) among ground states found:',len(pats))
