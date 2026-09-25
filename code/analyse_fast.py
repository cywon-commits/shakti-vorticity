import pickle, glob, math
import numpy as np
from vform2 import r4
rows=[]
for f in sorted(glob.glob('fast/f_*.pkl')):
    D=pickle.load(open(f,'rb')); n=D['n']; N=D['N']; x0,y0=D['x0'],D['y0']; cx,cy=D['center']
    E=D['emean']; V=D['vmean']
    # epsilon amplitude: max |<eps>| over the 3x3 vertices around the centre
    amp=max(abs(E[cx-x0+i,cy-y0+j]) for i in (-1,0,1) for j in (-1,0,1))
    ec=E[cx-x0,cy-y0]
    # v superstructure amplitude: Fourier components at (±pi/2, ±pi/2) over central disk of radius 0.25 n/sqrt2
    Rr=0.25*n/math.sqrt(2); acc={}; cnt=0
    for c in D['cells']:
        if math.hypot(c[0]+.5-cx,c[1]+.5-cy)<=Rr:
            cnt+=1; val=V[c[0]-x0,c[1]-y0]
            for q in [(1,1),(1,-1),(-1,1),(-1,-1),(1,0),(0,1),(2,0),(0,2),(2,2),(1,2),(2,1)]:
                acc[q]=acc.get(q,0)+val*np.exp(-1j*math.pi/2*(q[0]*c[0]+q[1]*c[1]))
    P={q:abs(a/cnt)**2 for q,a in acc.items()}
    mv=math.sqrt(P[(1,1)]+P[(1,-1)]+P[(-1,1)]+P[(-1,-1)])
    other=math.sqrt(max(P[q] for q in P if q not in [(1,1),(1,-1),(-1,1),(-1,-1)]))
    rows.append(dict(n=n,N=N,amp=amp,ec=ec,se=1/math.sqrt(N),var=D['hcenter'].var(),mv=mv,other=other,tot=int(D['tot'][0])))
rows.sort(key=lambda r:r['n'])
for r in rows: print('%3d N=%4d  eps_max=%.3f eps_c=%+.3f  Var h_c=%.2f  v-superstructure=%.3f  (largest other Fourier %.3f)  Σv=%+d'%(r['n'],r['N'],r['amp'],r['ec'],r['var'],r['mv'],r['other'],r['tot']))
pickle.dump(rows,open('fast/rows.pkl','wb'))
# fits
for lab,sel in [('all',lambda r:True),('n≡0',lambda r:r['n']%4==0),('n≡3',lambda r:r['n']%4==3),('n≡1',lambda r:r['n']%4==1),('n≡2',lambda r:r['n']%4==2)]:
    R=[r for r in rows if sel(r) and r['n']>=12]
    if len(R)<3: continue
    x=np.log([r['n'] for r in R])
    for key_ in ('amp','mv'):
        y=np.log([r[key_] for r in R]); s=np.polyfit(x,y,1)[0]
        print(lab,key_,'exponent %.3f'%(-s),end='  |  ')
    yv=[r['var'] for r in R]; print('dVar/dln n = %.3f'%np.polyfit(x,yv,1)[0])
