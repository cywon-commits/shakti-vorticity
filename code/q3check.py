import sys
sys.argv=['x','8','z3q3']
src=open('charge_pair.py').read().split("E0,_=solve_fix({})")[0]
src=src.replace("prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=600))","prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=600)); print(pulp.LpStatus[prob.status],end=' ')")
exec(src)
E0,_=solve_fix({}); print(E0)
vz=sorted([p for p in T.verts if len(T.inc[p])==3]); p1=vz[0]
for p2 in vz[1:8]:
    for s1 in (1,-1):
        E,sp=solve_fix({p1:('z3q3',s1),p2:('z3q3',-s1)}); print(p1,p2,s1,E)
