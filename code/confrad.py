import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, math
def conformal_radius(a, b=1.0, M=120):
    """conformal radius of rectangle [0,a]x[0,b] at its centre, by solving Laplace for u with u=ln|zeta-z0| on the boundary"""
    nx=int(round(a*M)); ny=int(round(b*M)); hx=a/nx; hy=b/ny
    z0=(a/2,b/2)
    N=(nx-1)*(ny-1); idx=lambda i,j:(i-1)*(ny-1)+(j-1)
    rows=[];cols=[];vals=[]; rhs=np.zeros(N)
    for i in range(1,nx):
        for j in range(1,ny):
            k=idx(i,j); rows.append(k); cols.append(k); vals.append(-2/hx**2-2/hy**2)
            for di,dj,hh in ((1,0,hx),(-1,0,hx),(0,1,hy),(0,-1,hy)):
                ii,jj=i+di,j+dj
                if 0<ii<nx and 0<jj<ny: rows.append(k); cols.append(idx(ii,jj)); vals.append(1/hh**2)
                else:
                    x,y=ii*hx,jj*hy; rhs[k]-=math.log(math.hypot(x-z0[0],y-z0[1]))/hh**2
    A=sp.csr_matrix((vals,(rows,cols)),shape=(N,N)); u=spl.spsolve(A,rhs)
    ic,jc=nx//2,ny//2
    return math.exp(u[idx(ic,jc)])
for a in (1,1.25,1.5,2,2.5,3,4,6,8):
    print(a, round(conformal_radius(a),4))
print('strip limit 2/pi =',round(2/math.pi,4))
