import io, contextlib, random
with contextlib.redirect_stdout(io.StringIO()):
    from placement2 import *
def rect_cells(w,h,sx,sy): return {(sx+i,sy+j) for i in range(w) for j in range(h)}
def house(w,sx,sy):
    c=rect_cells(w,w,sx,sy)
    # staircase roof: rows above with shrinking width (2 cells per row from each side), even widths kept
    k=0; ww=w
    while ww-4>=2:
        k+=1; ww-=4
    row=w; width=w
    while width>2:
        width-=2
        off=(w-width)//2
        c|={(sx+off+i,sy+row) for i in range(width)}
        row+=1
    return c
def annulus(w,sx,sy,hole=2):
    c=rect_cells(w,w,sx,sy); o=(w-hole)//2
    return c-{(sx+o+i,sy+o+j) for i in range(hole) for j in range(hole)}
def study(name,cells,cap=4000):
    C=Counter(); tots=Counter(); n=0
    for doms in enum_cells(cells):
        n+=1
        if n>cap: break
        I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
        C[('frustrated' if (bad or fr) else 'consistent')+f', {comps} comp']+=1
        if not (bad or fr):
            try:
                h=height_norm(doms,cells); f=[v_formula(d,h) for d in doms]; tots[abs(int(sum(f)))]+=1
            except Exception as e: tots['height multivalued']+=1
    print(f'{name:38s} cells={len(cells):3d} tilings checked={min(n,cap):5d} {dict(C)}  |Σv|: {dict(tots)}')
for sh in ((0,0),(1,0)):
    tag=f' shift{sh}'
    study('rectangle 4x6'+tag,rect_cells(4,6,*sh))
    study('rectangle 2x12'+tag,rect_cells(2,12,*sh))
    study('rectangle 6x8'+tag,rect_cells(6,8,*sh))
    study('house 6x6 + roof'+tag,house(6,*sh))
    study('house 8x8 + roof'+tag,house(8,*sh))
    study('annulus 6x6 minus 2x2 centre'+tag,annulus(6,*sh))
    study('annulus 8x8 minus 2x2 centre'+tag,annulus(8,*sh))
