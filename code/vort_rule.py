from vort_pairs import *
# defect (red) island of each tetromino = image of template edge ((0,1),(0,2))
def defect(A):
    return {e for e,(b,k) in A['tp'].items() if k=='per'} and None
T_DEF=((0,1),(0,2))
def info2(d):
    rs_,rl,s=tet_of(d); S=set(rs_)|set(rl)
    for f,refl in OPS:
        img={sq_img(f,x) for x in T_SQ}; imgS={sq_img(f,x) for x in T_SHORT}
        t=(min(x[0] for x in S)-min(x[0] for x in img), min(x[1] for x in S)-min(x[1] for x in img))
        if {(x[0]+t[0],x[1]+t[1]) for x in img}==S and {(x[0]+t[0],x[1]+t[1]) for x in imgS}==set(rs_):
            g=lambda p:(f(p)[0]+t[0],f(p)[1]+t[1]); break
    A=info(d); A['defect']=key((g(T_DEF[0]),g(T_DEF[1])))
    A['perim']={e for e,(b,k) in A['tp'].items() if k=='per'}
    return A
doms=[((i,j),(i+dx,j+dy)) for i in range(-2,2) for j in range(-2,2) for dx,dy in [(1,0),(0,1)]]
others=[((i,j),(i+dx,j+dy)) for i in range(-6,6) for j in range(-6,6) for dx,dy in [(1,0),(0,1)]]
C=Counter(); T=defaultdict(Counter)
for d in doms:
    A=info2(d)
    for o in others:
        if set(o)&set(d): continue
        B=info2(o); js,common,kinds=coupling(A,B)
        shared=A['perim']&B['perim']
        if not shared: continue
        rule={(-1)*(-1 if (e==A['defect'])!=(e==B['defect']) else 1) for e in shared}
        full=js
        C[(tuple(sorted(rule)),tuple(sorted(full)))]+=1
        for e in shared:
            T[('A-defect' if e==A['defect'] else 'A-normal')+' / '+('B-defect' if e==B['defect'] else 'B-normal')][A['tp'][e][0]*B['tp'][e][0]]+=1
print('rule vs full coupling:',dict(C))
print('per shared island:',{k:dict(v) for k,v in T.items()})
