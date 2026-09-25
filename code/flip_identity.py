"""Lemma 3 from the height formula: Sigma v of the flipped block is invariant (32 cases)."""
from vform2 import S_TAB, r4
def black(c): return (c[0]+c[1])%2==0
def step(p,q):
    dx,dy=q[0]-p[0],q[1]-p[1]
    L={(1,0):(p[0],p[1]),(-1,0):(q[0],q[1]-1),(0,1):(p[0]-1,p[1]),(0,-1):(q[0],q[1])}[(dx,dy)]
    return 1 if black(L) else -1
eps=lambda h:1 if (h//4)%2==0 else -1
bad=0
for cx in range(4):
    for cy in range(4):
        a=(cx-1,cy-1); b=(cx-1,cy); d=(cx,cy-1)
        for h0 in range(r4(a),r4(a)+8,4):
            vH=S_TAB[('H',a[0]%4,a[1]%4)]*eps(h0)+S_TAB[('H',b[0]%4,b[1]%4)]*eps(h0+step(a,b))
            vV=S_TAB[('V',a[0]%4,a[1]%4)]*eps(h0)+S_TAB[('V',d[0]%4,d[1]%4)]*eps(h0+step(a,d))
            bad+=vH!=vV
print('cases where a flip changes Sigma v:',bad,'of 32')
