"""Dipolar cost of an adjacent O pair relative to the two dominoes it replaces (point dipoles, cutoff 6)."""
import random, math, numpy as np, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import otet
from otet import tile_template, analyse, partner
from mono_tile import random_tiling
def spins_of(tiles,v):
    sp={}
    for k,t in enumerate(tiles):
        for e,(b,kind) in tile_template(t).items(): sp[e]=v[k]*b
    return sp
def energy(sp,centre,R=10,Rc=6.0,model='point',lam=0.7):
    E=[e for e in sp if math.hypot((e[0][0]+e[1][0])/2-centre[0],(e[0][1]+e[1][1])/2-centre[1])<R+Rc]
    pos=np.array([[(a[0]+b[0])/2,(a[1]+b[1])/2] for a,b in E]); d=np.array([[b[0]-a[0],b[1]-a[1]] for a,b in E],float); s=np.array([sp[e] for e in E],float)
    inner=np.hypot(pos[:,0]-centre[0],pos[:,1]-centre[1])<R
    tot=0.0
    for k in np.where(inner)[0]:
        dv=pos-pos[k]; r=np.hypot(dv[:,0],dv[:,1]); m=(r>0)&(r<=Rc)
        if model=='point':
            rv=dv[m]/r[m,None]; J=(d[m]@d[k]-3*(rv@d[k])*np.sum(d[m]*rv,axis=1))/r[m]**3
        else:
            J=np.zeros(m.sum()); P=pos[m]; Dm=d[m]
            for sa in (1,-1):
                for sb in (1,-1):
                    pa=pos[k]+sa*lam/2*d[k]; pb=P+sb*lam/2*Dm; J+=sa*sb/np.hypot(*(pa-pb).T)
            J/=lam*lam
        w=np.where(inner[m],0.5,1.0)     # pairs inside counted once
        tot+=s[k]*np.sum(w*J*s[m])
    return tot
from core import rect
rng=random.Random(4); L=18; board={(i,j) for i in range(L) for j in range(L)}
for model in ('point','dumbbell'):
    dEs=[]
    for trial in range(40):
        c=(L//2+rng.randrange(-2,2),L//2+rng.randrange(-2,2)); c=c if (c[0]+c[1])%2==0 else (c[0]+1,c[1])
        pc=partner(c)
        # find c2 adjacent to c such that {c,c2} and {pc,p(c2)} are dominoes
        adj=lambda a,b: abs(a[0]-b[0])+abs(a[1]-b[1])==1
        c2=next((x for x in [(c[0]+2,c[1]+1),(c[0]+1,c[1]+2),(c[0]-1,c[1]+2),(c[0]+1,c[1]),(c[0],c[1]+1)] if (x[0]+x[1])%2==1 and partner(x) in board and not {x,partner(x)}&{c,pc}),None)
        if c2 is None: continue
        p2=partner(c2); Ocells={c,pc,c2,p2}
        rest=random_tiling(board-Ocells,rng)
        B=rest+[('O',)+tuple(sorted((c,pc))),('O',)+tuple(sorted((c2,p2)))]
        touch=[t for t in rest if any(adj(x,y) for x in t[1:] for y in Ocells)]
        W=set(Ocells)|{x for t in touch for x in t[1:]}
        outside=[t for t in rest if t not in touch]
        _p=open('placement2.py').read(); exec(_p[_p.index('def enum_cells'):_p.index('for n in range(1,6)')],globals())
        alts=[outside+[('D',)+tuple(sorted(d)) for d in dd] for dd in enum_cells(W)]
        if not alts: continue
        cb,fb,nb_,vb=analyse(B)
        if cb or fb: continue
        cen=(np.mean([(p[0]+.5) for x in Ocells for p in rect(x)]),np.mean([(p[1]+.5) for x in Ocells for p in rect(x)]))
        EAs=[]
        for A in alts:
            ca,fa,na,va=analyse(A)
            if ca or fa: continue
            far=max(outside,key=lambda t:abs(t[1][0]-c[0])+abs(t[1][1]-c[1])); ka=A.index(far); kb=B.index(far); g=va[ka]*vb[kb]; vbb={k:g*x for k,x in vb.items()}
            EAs.append((energy(spins_of(A,va),cen,model=model),vbb))
        if not EAs: continue
        EA=min(e for e,_ in EAs); vbb=EAs[0][1]
        dEs.append(energy(spins_of(B,vbb),cen,model=model)-EA)
        continue
        # align global sign on the first common (far) tile
        g=va[0]*vb[0]; vb={k:g*x for k,x in vb.items()}
        cen=np.mean([[(p[0]+.5) for p in rect(x)] for x in block]),np.mean([[(p[1]+.5) for p in rect(x)] for x in block])
        dEs.append(energy(spins_of(B,vb),cen,model=model)-energy(spins_of(A,va),cen,model=model))
    dEs=np.array(dEs); print(model,'adjacent O pair: dE = %.2f +- %.2f D (min %.2f, max %.2f, n=%d)'%(dEs.mean(),dEs.std(),dEs.min(),dEs.max(),len(dEs)))
