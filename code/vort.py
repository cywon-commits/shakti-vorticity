"""Vorticity rule for L-tetrominoes on the Shakti lattice, read off the user's figure.
v = +1 : counter-clockwise circulation (orange), v = -1 : clockwise (blue). Time reversal = flip all arrows."""
from core import *
from collections import Counter, defaultdict
# ---- template: orange L (v=+1). squares (0,0),(1,0) = short horizontal arm ; (0,1),(0,2) = long vertical arm
T_PER=[((0,0),(1,0)),((1,0),(2,0)),((2,0),(2,1)),((2,1),(1,1)),((1,1),(1,2)),((1,2),(1,3)),((1,3),(0,3)),((0,3),(0,2)),
       ((0,1),(0,2)),   # red: against the circulation
       ((0,1),(0,0))]
T_INT=[((1,1),(0,1))]
T_EXT=[((0,3),(0,4)),((-1,3),(0,3)),((0,2),(-1,2)),((-1,1),(0,1)),((2,2),(2,1)),((2,1),(3,1))]
T_FREE=[((0,0),(0,-1)),((2,0),(2,-1)),((1,3),(2,3))]
T_SQ={(0,0),(1,0),(0,1),(0,2)}; T_SHORT={(0,0),(1,0)}
# ---- D4 on vertex coordinates (squares transform via their centre)
OPS=[]
for refl in (False,True):
    for r in range(4):
        def f(p,r=r,refl=refl):
            x,y=p
            if refl: x=-x
            for _ in range(r): x,y=-y,x
            return (x,y)
        OPS.append((f,refl))
def sq_img(f,s):
    cx,cy=f((s[0]+.5,s[1]+.5)); return (int(math.floor(cx)),int(math.floor(cy)))
def key(e):
    a,b=e; return (a,b) if a<b else (b,a)
def sgn(e):  # +1 if directed along canonical order
    a,b=e; return 1 if a<b else -1
REMOVED=set()
def removed_edges_near(sqs):
    pass
def is_island(e):
    """edge exists in Shakti lattice iff it is not the middle edge of a plaquette"""
    (x1,y1),(x2,y2)=key(e)
    if x1==x2:   # vertical edge x=x1, y1..y1+1 : middle of horizontal plaquette {(x1-1,y1),(x1,y1)}?
        a,b=math.floor((x1-1)/2),math.floor(y1/2)
        return not ((a+b)%2==0 and x1%2==1)
    else:        # horizontal edge y=y1 : middle of vertical plaquette {(x1,y1-1),(x1,y1)}
        a,b=math.floor(x1/2),math.floor((y1-1)/2)
        return not ((a+b)%2==1 and y1%2==1)
def template_for(short_sq, long_sq):
    """return dict canonical-edge -> base sign (for v=+1) and set of free edges, for a tetromino in the lattice"""
    S=set(short_sq)|set(long_sq)
    for f,refl in OPS:
        img={sq_img(f,s) for s in T_SQ}; imgS={sq_img(f,s) for s in T_SHORT}
        # translation
        t=(min(s[0] for s in S)-min(s[0] for s in img), min(s[1] for s in S)-min(s[1] for s in img))
        if {(s[0]+t[0],s[1]+t[1]) for s in img}==S and {(s[0]+t[0],s[1]+t[1]) for s in imgS}==set(short_sq):
            g=lambda p:(f(p)[0]+t[0],f(p)[1]+t[1])
            rev=-1 if refl else 1
            out={}
            for grp,kind in ((T_PER,'per'),(T_INT,'int'),(T_EXT,'ext')):
                for a,b in grp:
                    e=(g(a),g(b)); assert is_island(e),('missing island',e,kind)
                    out[key(e)]=(sgn(e)*rev,kind)
            free={key((g(a),g(b))) for a,b in T_FREE}
            for e in free: assert is_island(e)
            return out,free
    raise ValueError('no template match')
def tet_of(d):
    s=short_cell(*d); l=d[0] if s==d[1] else d[1]
    return rect(s), rect(l), s
if __name__=='__main__':
    # sanity: every orientation of L occurring in the lattice matches the template and islands exist
    import random
    cnt=Counter()
    for i in range(-6,6):
        for j in range(-6,6):
            for dx,dy in [(1,0),(0,1)]:
                d=((i,j),(i+dx,j+dy)); rs_,rl,s=tet_of(d); tp,fr=template_for(rs_,rl); cnt[len(tp)]+=1
    print('all L placements matched; edges per template',cnt)
