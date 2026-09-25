import pickle, glob, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from core import rect, pcen
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
def amp(D):
    L=D['L']; sx,sy=D['shift']; V=D['vmean']; acc={}; cnt=0; c0=L/2
    for i in range(L):
        for j in range(L):
            if math.hypot(i+.5-c0,j+.5-c0)<=L/4:
                cnt+=1; x,y=sx+i,sy+j
                for q in [(1,1),(1,-1),(-1,1),(-1,-1),(1,0),(0,1),(2,0),(0,2),(2,2),(1,2),(2,1)]:
                    acc[q]=acc.get(q,0)+V[i,j]*np.exp(-1j*math.pi/2*(q[0]*x+q[1]*y))
    P={q:abs(a/cnt)**2 for q,a in acc.items()}
    main=math.sqrt(sum(P[q] for q in [(1,1),(1,-1),(-1,1),(-1,-1)]))
    other=math.sqrt(max(P[q] for q in P if q not in [(1,1),(1,-1),(-1,1),(-1,-1)]))
    return main,other
rows=[]
for f in sorted(glob.glob('sq/sq_*.pkl')):
    D=pickle.load(open(f,'rb')); m,o=amp(D); L=D['L']
    rows.append(dict(L=L,cls='black' if sum(D['shift'])%2==0 else 'white',amp=m,other=o,var=D['hcenter'].var(),ec=D['emean'][L//2,L//2],tot=sorted(set(D['tot'].tolist())),N=D['N']))
rows.sort(key=lambda r:(r['cls'],r['L']))
for r in rows: print('%s L=%2d N=%d amp=%.3f other=%.3f Var=%.2f eps_c=%+.3f Σv=%s'%(r['cls'],r['L'],r['N'],r['amp'],r['other'],r['var'],r['ec'],r['tot']))
for cls in ('black','white'):
    for par in (0,2):
        R=[r for r in rows if r['cls']==cls and r['L']%4==par and r['L']>=12]
        if len(R)>2:
            s=np.polyfit(np.log([r['L'] for r in R]),np.log([r['amp'] for r in R]),1)[0]; print(cls,'L≡%d'%par,'exponent %.3f'%(-s))
    R=[r for r in rows if r['cls']==cls and r['L']>=12]
    print(cls,'all exponent %.3f'%(-np.polyfit(np.log([r['L'] for r in R]),np.log([r['amp'] for r in R]),1)[0]))
pickle.dump(rows,open('sq/rows.pkl','wb'))
