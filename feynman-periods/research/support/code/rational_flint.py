"""Small exact rational-function wrapper around local FLINT multivariate polynomials."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import fmpq,fmpq_mpoly_ctx
CTX=fmpq_mpoly_ctx.get(('m12','m22','m32','m42'))
PZERO=CTX.constant(0);PONE=CTX.constant(1)
def po(p):return [[*map(int,e),str(c)] for e,c in sorted(p.to_dict().items())]
def pi(data):return CTX.from_dict({tuple(v[:4]):fmpq(v[4]) for v in data})
class RF:
    __slots__=('n','d')
    def __init__(self,n=0,d=1):
        if isinstance(n,RF) and d==1:self.n,self.d=n.n,n.d;return
        n=CTX.constant(n) if isinstance(n,(int,fmpq)) else n
        d=CTX.constant(d) if isinstance(d,(int,fmpq)) else d
        assert d!=0
        if not n:self.n,self.d=PZERO,PONE;return
        g=n.gcd(d);n=n/g;d=d/g;c=d.coeffs()[0]
        self.n,self.d=n/c,d/c
    def __add__(self,b):
        b=RF(b)
        if not self.n:return b
        if not b.n:return self
        if self.d==b.d:return RF(self.n+b.n,self.d)
        g=self.d.gcd(b.d);a=self.d/g;c=b.d/g
        return RF(self.n*c+b.n*a,a*b.d)
    __radd__=__add__
    def __neg__(self):return RF(-self.n,self.d)
    def __sub__(self,b):return self+-RF(b)
    def __rsub__(self,b):return RF(b)+-self
    def __mul__(self,b):
        b=RF(b)
        if not self.n or not b.n:return RF(0)
        g=self.n.gcd(b.d);h=b.n.gcd(self.d)
        return RF((self.n/g)*(b.n/h),(self.d/h)*(b.d/g))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=RF(b);assert b.n!=0
        return self*RF(b.d,b.n)
    def __rtruediv__(self,b):return RF(b)/self
    def __eq__(self,b):
        b=RF(b);return self.n*b.d==b.n*self.d
    def __bool__(self):return bool(self.n)
    def diff(self,k):return RF(self.n.derivative(k)*self.d-self.n*self.d.derivative(k),self.d**2)
    def theta(self,k):return RF(CTX.gen(k))*self.diff(k)
    def encode(self):return {'n':po(self.n),'d':po(self.d)}
    @classmethod
    def decode(cls,v):return cls(pi(v['n']),pi(v['d']))
    def __repr__(self):return '('+str(self.n)+')/('+str(self.d)+')'
X=[RF(x) for x in CTX.gens()]
def read_matrix(path):
    return [[RF.decode(v) for v in row] for row in json.loads(Path(path).read_text())['matrix']]
def save(path,obj):
    p=Path(path);p.parent.mkdir(exist_ok=True);q=p.with_suffix('.tmp')
    q.write_text(json.dumps(obj,ensure_ascii=False,separators=(',',':'))+'\n',encoding='utf-8');q.replace(p)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),RF(0))
