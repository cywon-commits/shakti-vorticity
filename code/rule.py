exec(open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0].split("def pale")[0])
from collections import Counter
C=Counter()
for i in range(-6,6):
  for j in range(-6,6):
    if not black((i,j)): continue
    blk=inv[(i,j)][:2]
    part=next(c for c in [(i+1,j+1),(i-1,j-1),(i+1,j-1),(i-1,j+1)] if c in inv and inv[c][:2]==blk)
    for t,(dx,dy) in {'E':(1,0),'W':(-1,0),'N':(0,1),'S':(0,-1)}.items():
        w=(i+dx,j+dy); d=((i,j),w); sc=short_cell(d)
        wb=inv[w][:2]; wpart=next(c for c in [(w[0]+1,w[1]+1),(w[0]-1,w[1]-1),(w[0]+1,w[1]-1),(w[0]-1,w[1]+1)] if c in inv and inv[c][:2]==wb)
        C[(t,(part[0]-i,part[1]-j),(wpart[0]-w[0],wpart[1]-w[1]),'short=B' if sc==d[0] else 'short=W')]+=1
for k in sorted(C): print(k,C[k])
bad=0;tot=0
from collections import defaultdict
D=defaultdict(Counter)
for i in range(-6,6):
  for j in range(-6,6):
    if not black((i,j)): continue
    for t,(dx,dy) in {'E':(1,0),'W':(-1,0),'N':(0,1),'S':(0,-1)}.items():
        w=(i+dx,j+dy); tot+=1; bad+= not isL(rs[inv[(i,j)]],rs[inv[w]])
        # geometric relation: does W plaquette touch B plaquette end-on (at B's short end)?
        rb=sorted(rs[inv[(i,j)]]); sc=short_cell(((i,j),w))
        D[(t,(i%2,j%2) if False else ((i+j)//2)%2, i%2)]['B' if sc==(i,j) else 'W']+=1
print(bad,tot)
for k in sorted(D): print(k,dict(D[k]))
