"""Open Shakti sample cut as an axis-aligned box: exact vertex-model ground state (ILP), unhappy vertices -> dimers."""
import sys, itertools, io, contextlib, pulp
from collections import defaultdict, Counter
with contextlib.redirect_stdout(io.StringIO()):
    from vort import is_island
from shakti_torus import classify
from axis_region import CELLS, shakti_rect
from exact_v import aztec_cells
COST={'z1':0,'z2ok':0,'z3happy':0,'z3unhappy':1,'z4I':0,'z4II':3,'z4q2':3,'z3q3':3,'z2q':3,'z4q4':20}
def solve(x0,y0,X,Y,force=None):
    """force: plaquette cells that must each carry an unhappy vertex (all plaquettes in tetrominoes)"""
    isl=[]
    for x in range(x0,x0+X+1):
        for y in range(y0,y0+Y+1):
            for e in (((x,y),(x+1,y)),((x,y),(x,y+1))):
                if e[1][0]<=x0+X and e[1][1]<=y0+Y and is_island(e): isl.append(e)
    inc=defaultdict(list)
    for k,(a,b) in enumerate(isl): inc[a].append((k,'start',a,b)); inc[b].append((k,'end',a,b))
    prob=pulp.LpProblem('open',pulp.LpMinimize); sv={k:pulp.LpVariable(f's{k}',cat='Binary') for k in range(len(isl))}
    obj=[]; Y_={}
    def lab(L,bits):
        ins=[];dirs=[]
        for (k,role,a,b),bt in zip(L,bits):
            ins.append((bt==1 and role=='end') or (bt==-1 and role=='start')); dirs.append((b[0]-a[0],b[1]-a[1]))
        z=len(L); q=sum(1 if i else -1 for i in ins)
        return 'z1' if z==1 else classify(z,q,ins,dirs)
    for p,L in inc.items():
        rows=[(bits,lab(L,bits)) for bits in itertools.product((1,-1),repeat=len(L))]
        ys=[pulp.LpVariable(f'y{p[0]}_{p[1]}_{a}'.replace('-','m'),cat='Binary') for a in range(len(rows))]
        prob+=pulp.lpSum(ys)==1; obj+= [COST[l]*v for (b,l),v in zip(rows,ys)]
        for j,(k,_,_,_) in enumerate(L): prob+= sv[k]==pulp.lpSum(v for (b,l),v in zip(rows,ys) if b[j]==1)
        Y_[p]=(rows,ys)
    if force:
        unh={p:pulp.lpSum(v for (b,l),v in zip(*Y_[p]) if l=='z3unhappy') for p in Y_}
        for c in force:
            r=CELLS[c]; xs=[q[0] for q in r]; ys=[q[1] for q in r]; X0,X1,Y0,Y1=min(xs),max(xs)+1,min(ys),max(ys)+1
            bd=[(x,y) for x in range(X0,X1+1) for y in range(Y0,Y1+1) if (x in (X0,X1) or y in (Y0,Y1)) and (x,y) in unh]
            prob+= pulp.lpSum(unh[p] for p in bd)>=1
        prob+= pulp.lpSum(unh.values())>=len(force)//2
    prob+=pulp.lpSum(obj); prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=600))
    spins={k:(1 if sv[k].value()>.5 else -1) for k in sv}
    labs={p:lab(inc[p],[spins[k] for k,_,_,_ in inc[p]]) for p in inc}
    solve.last=dict(isl=isl,spins=spins)
    return pulp.value(prob.objective), labs
if __name__=='__main__':
    n=int(sys.argv[1])
    for pl in ('II','I'):
        A,_=aztec_cells(n,pl); pts=[p for c in A for p in CELLS[c]]
        x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
        box=shakti_rect(x0,y0,x1-x0,y1-y0)
        E,labs=solve(x0,y0,x1-x0,y1-y0)
        un=[p for p,l in labs.items() if l=='z3unhappy']
        # plaquettes touched by each unhappy vertex
        def corners(c):
            r=CELLS[c]; xs=[q[0] for q in r]; ys=[q[1] for q in r]
            return {(x,y) for x in range(min(xs),max(xs)+2) for y in range(min(ys),max(ys)+2)}
        cov=Counter(); dimers=[]
        for p in un:
            pl_=[c for c in box if p in corners(c) and (p[0],p[1]) not in {(q[0]+.5,q[1]+.5) for q in CELLS[c]}]
            dimers.append(pl_)
            for c in pl_: cov[c]+=1
        bad=Counter(l for l in labs.values() if COST[l]>1)
        print(f'n={n} {pl}: box {x1-x0}x{y1-y0}, complete plaquettes {len(box)} (Aztec {len(A)}), energy {E}, unhappy {len(un)}, other costly {dict(bad)}')
        print('   plaquettes with >=1 unhappy vertex:',sum(1 for c in box if cov[c]>0),' of which in Aztec:',sum(1 for c in A if cov[c]>0),' edge extras with unhappy:',sum(1 for c in box-A if cov[c]>0))
