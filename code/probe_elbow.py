# which vertex of an L template is unhappy? and are all template vertices otherwise OK?
from shakti_torus import *
T=Torus(8)
c=T.cands[0]; t=T.tmpl(c[2],c[3])
# vertices whose all islands are in template
for p in T.verts:
    ks=[k for k,_,_,_ in T.inc[p]]
    if all(k in t for k in ks):
        z,q,ins,dirs=T.vstate(p,t); print(p,classify(z,q,ins,dirs))
print('squares',sorted(c[2]|c[3]),'elbow',elbow(c[2]|c[3]))
