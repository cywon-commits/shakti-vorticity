"""Free minimisation of the Shakti vertex model on a torus with a +/- pair of charge defects of a given type forced."""
import sys, math, time
from mono_ilp import *
P=int(sys.argv[1]); typ=sys.argv[2]; cost=dict(COST0)
for a in sys.argv[3:]:
    k,v=a.split('='); cost[k]=float(v)
M=Model(P); T=M.T
def solve_fix(fix):
    prob=pulp.LpProblem('s',pulp.LpMinimize); y={}; obj=[]
    for p,rows in M.tab.items():
        for a,(bits,lab) in enumerate(rows):
            y[p,a]=pulp.LpVariable(f'y_{p[0]}_{p[1]}_{a}',cat='Binary'); obj.append(cost[lab]*y[p,a])
        allowed=range(len(rows))
        if p in fix:
            lab0,qs=fix[p]
            allowed=[a for a,(bits,l) in enumerate(rows) if l==lab0 and qs*T.vstate(p,{k:b for (k,_,_,_),b in zip(T.inc[p],bits)})[1]>0]
        prob+= pulp.lpSum(y[p,a] for a in allowed)==1
    sv={k:pulp.LpVariable(f's_{k}',cat='Binary') for k in range(len(T.isl))}
    for p,rows in M.tab.items():
        for j,(k,_,_,_) in enumerate(T.inc[p]):
            prob+= sv[k]==pulp.lpSum(y[p,a] for a,(bits,l) in enumerate(rows) if bits[j]==1)
    prob+= pulp.lpSum(obj); prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=600))
    spins={k:(1 if sv[k].value()>.5 else -1) for k in sv}
    return pulp.value(prob.objective), spins
E0,_=solve_fix({}); print('ground state E',E0,'(plaquettes/2 =',len(T.plaqs)/2,')',flush=True)
z={'z4q2':4,'z3q3':3,'z2q':2,'z4II':4}[typ]
vz=sorted([p for p in T.verts if len(T.inc[p])==z])
p1=vz[0]
def td(a,b):
    dx=abs(a[0]-b[0]); dy=abs(a[1]-b[1]); return math.hypot(min(dx,P-dx),min(dy,P-dy))
out=[]; seen=set()
for p2 in sorted(vz[1:],key=lambda q:td(p1,q)):
    r=round(td(p1,p2),2)
    key_=(r,)
    if key_ in seen: continue
    seen.add(key_)
    best=None
    for s1 in ((1,-1) if typ!='z4II' else (1,)):
        fix={p1:(typ,s1 if typ!='z4II' else 1),p2:(typ,-s1 if typ!='z4II' else 1)} if typ!='z4II' else None
        if typ=='z4II':
            # neutral: fix label only
            fix={p1:(typ,0),p2:(typ,0)}
        E,sp=solve_fix(fix)
        if best is None or E<best[0]: best=(E,sp,s1)
    E,sp,s1=best
    labs=Counter(classify(*T.vstate(p,sp)) for p in T.verts)
    print(f'{typ} pair r={r:5.2f}: excess E = {E-E0:5.1f}  vertices {dict(labs)}',flush=True)
    out.append((r,E-E0,p1,p2,s1,sp))
    import pickle; pickle.dump(dict(P=P,typ=typ,cost=cost,E0=E0,out=out),open(f'mono/pair{P}_{typ}.pkl','wb'))
    if r>P/2: break
