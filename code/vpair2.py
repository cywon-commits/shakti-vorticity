import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vort_rule import *
C=Counter(); cls={}
for x in range(-4,4):
    for y in range(-4,4):
        J={}
        for fill in ('h','v'):
            if fill=='h': d1,d2=((x,y),(x+1,y)),((x,y+1),(x+1,y+1))
            else: d1,d2=((x,y),(x,y+1)),((x+1,y),(x+1,y+1))
            J[fill]=coupling(info2(d1),info2(d2))[0].pop()
        key_=((x+y)%4,(x-y)%4)
        C[(J['h'],J['v'])]+=1; cls.setdefault(key_,set()).add((J['h'],J['v']))
print(C); print(cls)
