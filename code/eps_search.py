from fractions import Fraction as F
import itertools
from eps_formula import vals
r={m:vals[m+1]/vals[m] for m in range(1,8)}
d={m:r[m+1]/r[m] for m in range(1,7)}
# model: d_m = K * prod_j ((4m+4+j)/(4m+j))^{c_j} * prod_j ((2m+2+j)/(2m+j))^{e_j}, j=1..3 (odd for 2m: j=1)
best=[]
rng=range(-6,7)
for c in itertools.product(rng,repeat=4):   # j=0..3 for 4m+j ; (4m+j) with j=0 -> m factors
    def P(m):
        p=F(1)
        for j,cj in enumerate(c):
            if cj: p*=F(4*m+4+j,4*m+j)**cj if 4*m+j>0 else 1
        return p
    K=d[1]/P(1)
    if all(d[m]==K*P(m) for m in range(2,7)): best.append((c,K))
print('solutions with 4m+j factors:',best[:5])
