"""Phase map: exact arg Psi over the whole region vs pi<h>/4 from the exact mean height (limit-shape prediction)."""
import sys, pickle, math, cmath, numpy as np
from exact_psi import psi_map
from exact_meanh import mean_height
from exact_v import *
out={}
for name,cells in [('Aztec I n=24',aztec_cells(24,'I')[0]),('Aztec II n=24',aztec_cells(24,'II')[0]),('square L=40',{(i,j) for i in range(40) for j in range(40)}),
                   ('rectangle 56x28',{(i,j) for i in range(56) for j in range(28)})]:
    R=Region(cells); s,hs=start_vertex(cells)
    P,E=psi_map(R,s,hs,sorted(cells)); H=mean_height(R)
    hc={c:np.mean([H[(c[0]+a,c[1]+b)] for a in (0,1) for b in (0,1)]) for c in cells}
    d=[((cmath.phase(P[c])-math.pi*hc[c]/4+math.pi)%(2*math.pi))-math.pi for c in P if abs(P[c])>0.05]
    out[name]=dict(cells=cells,P=P,hc=hc)
    print(name,'cells with |Psi|>0.05:',len(d),' max |arg Psi - pi<h>/4| = %.3f rad, rms %.3f'%(max(abs(x) for x in d),np.sqrt(np.mean(np.square(d)))),flush=True)
pickle.dump(out,open('exact/phase_map.pkl','wb'))
