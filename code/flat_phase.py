import sys, math, cmath, numpy as np
from vform2 import S_TAB, r4
black=lambda c:(c[0]+c[1])%2==0
def partner_fn(al,be,par):
    def P(c):
        x,y=c; a=x-((x-al)%2); b=y-((y-be)%2)          # block lower-left cell
        horiz=(((a-al)//2+(b-be)//2)%2==0)!=par
        if horiz: return (a+(1-(x-a)), y)               # pair along x within the block row
        return (x, b+(1-(y-b)))
    return P
def heights(P,W):
    # BFS over vertices in [-W,W]^2 from origin, h(0,0)=r4(0,0)
    from collections import deque
    h={(0,0):r4((0,0))}; Q=deque([(0,0)])
    while Q:
        u=Q.popleft()
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            v=(u[0]+d[0],u[1]+d[1])
            if v in h or max(abs(v[0]),abs(v[1]))>W: continue
            dx,dy=d
            if (dx,dy)==(1,0): L,R=(u[0],u[1]),(u[0],u[1]-1)
            elif (dx,dy)==(-1,0): L,R=(v[0],v[1]-1),(v[0],v[1])
            elif (dx,dy)==(0,1): L,R=(u[0]-1,u[1]),(u[0],u[1])
            else: L,R=(v[0],v[1]),(v[0]-1,v[1])
            s=1 if black(L) else -1
            if P(L)==R: s=s-4 if s>0 else s+4
            h[v]=h[u]+s; Q.append(v)
    return h
def analyse(al,be,par,off=0,W=12):
    P=partner_fn(al,be,par); h={k:v+off for k,v in heights(P,W).items()}
    eps=lambda p:1 if (h[p]//4)%2==0 else -1
    cells=[(x,y) for x in range(-8,8) for y in range(-8,8)]
    psi=np.mean([0.5*sum(cmath.exp(1j*math.pi*r4((x+a,y+b))/4)*eps((x+a,y+b)) for a in (0,1) for b in (0,1)) for (x,y) in cells])
    hm=np.mean([h[(x,y)]-r4((x,y)) for x in range(-8,9) for y in range(-8,9)])   # mean offset above residue
    hbar=np.mean([h[(x,y)] for x in range(-8,9) for y in range(-8,9)])
    return psi,hbar,P,h
if __name__=='__main__':
    for par in (False,True):
        for al in (0,1):
            for be in (0,1):
                for off in (0,4):
                    psi,hb,_,_=analyse(al,be,par,off)
                    print(f'par={int(par)} (al,be)=({al},{be}) off={off}: |Psi|={abs(psi):.4f} argPsi/(pi/4)={cmath.phase(psi)/(math.pi/4):+.3f}  <h>={hb:.3f}  pi<h>/4 /(pi/4) mod 8={hb%8:.3f}')
