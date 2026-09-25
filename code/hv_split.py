"""Split of the central <v> order between horizontal (black-cell) and vertical (white-cell) plaquettes."""
import pickle, glob, math
import numpy as np
def split_aztec(D):
    n=D['n']; N=D['N']; x0,y0=D['x0'],D['y0']; cx,cy=D['center']; V=D['vmean']; R=0.25*n/math.sqrt(2)
    s={0:[],1:[]}
    for c in D['cells']:
        if math.hypot(c[0]+.5-cx,c[1]+.5-cy)<=R: s[(c[0]+c[1])%2].append(V[c[0]-x0,c[1]-y0])
    return s,N
def split_square(D):
    L=D['L']; N=D['N']; sx,sy=D['shift']; V=D['vmean']; s={0:[],1:[]}
    for i in range(L):
        for j in range(L):
            if math.hypot(i+.5-L/2,j+.5-L/2)<=L/4: s[(sx+i+sy+j)%2].append(V[i,j])
    return s,N
def amps(s,N,neff=None):
    out={}
    for k in (0,1):
        a=np.array(s[k]); m2=np.mean(a**2)-(1.0/(neff or N))   # subtract sampling-noise bias
        out[k]=math.sqrt(max(m2,0))
    H,Vv=out[0],out[1]
    return H,Vv,(H**2-Vv**2)/(H**2+Vv**2)
rows=[]
for d,lab in (('fast','Aztec I'),('fast2','Aztec II')):
    for f in glob.glob(d+'/f_*.pkl'):
        D=pickle.load(open(f,'rb'))
        if D['N']<600: continue
        s,N=split_aztec(D); H,Vv,asym=amps(s,N); rows.append(dict(geom=lab,size=D['n'],H=H,V=Vv,asym=asym))
for f in glob.glob('sq/sq_*.pkl'):
    D=pickle.load(open(f,'rb')); s,N=split_square(D); H,Vv,asym=amps(s,N,neff=N/10)
    rows.append(dict(geom='Square, corner '+('black' if sum(D['shift'])%2==0 else 'white'),size=D['L'],H=H,V=Vv,asym=asym))
rows.sort(key=lambda r:(r['geom'],r['size']))
pickle.dump(rows,open('hv_rows.pkl','wb'))
for g in sorted({r['geom'] for r in rows}):
    print('==',g)
    for r in rows:
        if r['geom']==g: print('  size %2d (mod4=%d)  A_H=%.3f  A_V=%.3f  asym=%+.2f'%(r['size'],r['size']%4,r['H'],r['V'],r['asym']))
