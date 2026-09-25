from vort import *
def plaq_edges(r):
    E=set()
    for (x,y) in r:
        for e in [((x,y),(x+1,y)),((x,y+1),(x+1,y+1)),((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]: E.add(key(e))
    return E
def orient(r): return 'H' if len({p[1] for p in r})==1 else 'V'
def info(d):
    rs_,rl,s=tet_of(d); tp,fr=template_for(rs_,rl)
    return dict(S=rs_,L=rl,tp=tp,free=fr,cells=set(d))
def coupling(A,B):
    common=set(A['tp'])&set(B['tp'])
    js={A['tp'][e][0]*B['tp'][e][0] for e in common}
    kinds=tuple(sorted((A['tp'][e][1],B['tp'][e][1]) for e in common))
    return js,common,kinds
def contact(A,B):
    out=[]
    for na,ra in (('s',A['S']),('l',A['L'])):
        for nb,rb in (('s',B['S']),('l',B['L'])):
            sh=plaq_edges(ra)&plaq_edges(rb)
            if sh: out.append(na+nb+('∥' if orient(ra)==orient(rb) else '⊥'))
    return tuple(sorted(out))
if __name__=='__main__':
    doms=[((i,j),(i+dx,j+dy)) for i in range(-2,2) for j in range(-2,2) for dx,dy in [(1,0),(0,1)]]
    others=[((i,j),(i+dx,j+dy)) for i in range(-6,6) for j in range(-6,6) for dx,dy in [(1,0),(0,1)]]
    T=defaultdict(Counter)
    for d in doms:
        A=info(d)
        for o in others:
            if set(o)&set(d): continue
            B=info(o); js,common,kinds=coupling(A,B)
            ct=contact(A,B)
            if not common and not ct: continue
            res='conflict' if len(js)>1 else ({1:'same',-1:'opposite'}[js.pop()] if js else 'no shared island')
            T[ct if ct else ('no plaquette contact',)][res]+=1
    for k,v in sorted(T.items(),key=lambda x:str(x[0])): print(k,dict(v))
