import pickle, math, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
def cr_square(z, M=96):
    n=M; h=1.0/n; N=(n-1)**2; idx=lambda i,j:(i-1)*(n-1)+(j-1)
    rows=[];cols=[];vals=[]; rhs=np.zeros(N)
    for i in range(1,n):
        for j in range(1,n):
            k=idx(i,j); rows.append(k); cols.append(k); vals.append(-4)
            for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
                ii,jj=i+di,j+dj
                if 0<ii<n and 0<jj<n: rows.append(k); cols.append(idx(ii,jj)); vals.append(1)
                else: rhs[k]-=math.log(math.hypot(ii*h-z[0],jj*h-z[1]))
    u=spl.spsolve(sp.csr_matrix((vals,(rows,cols)),shape=(N,N)),rhs)
    return math.exp(u[idx(int(round(z[0]*n)),int(round(z[1]*n)))])
out=[]
for sh in ('00','10'):
    D=pickle.load(open(f'sq/sq_48_{sh}.pkl','rb')); V=D['vmean']; L=48; sx,sy=D['shift']
    for t in (0.5,0.4,0.3,0.22,0.16):   # window centre along the horizontal midline at x = t*L
        cx=t*L; cy=L/2; acc={}; cnt=0
        for i in range(L):
            for j in range(L):
                if abs(i+.5-cx)<=4 and abs(j+.5-cy)<=4:
                    cnt+=1
                    for q in [(1,1),(1,-1),(-1,1),(-1,-1)]:
                        acc[q]=acc.get(q,0)+V[i,j]*np.exp(-1j*math.pi/2*(q[0]*(sx+i)+q[1]*(sy+j)))
        A=math.sqrt(sum(abs(a/cnt)**2 for a in acc.values()))
        out.append((sh,t,A,cr_square((t,0.5))))
for sh in ('00','10'):
    R=[o for o in out if o[0]==sh]; A0,r0=R[0][2],R[0][3]
    print('corner','black' if sh=='00' else 'white')
    for _,t,A,r in R: print('  x/L=%.2f  amplitude %.3f  ratio %.3f   predicted (rho/rho_c)^(-1/4) %.3f'%(t,A,A/A0,(r/r0)**-0.25))
