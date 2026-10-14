"""Exact four-variable arithmetic; restricted AST parsing of published data."""
from pathlib import Path
import ast, hashlib, json, re
from sympy import QQ
from sympy.polys.fields import field

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'external/bonn-banana-2025'
OUT = ROOT / 'data/banana_four_parameter'
FF, *XX = field('m12,m22,m32,m42', QQ)
RR = FF.ring
NAMES = dict(zip(('m12','m22','m32','m42'), XX))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def split_list(text):
    text = text.strip()
    assert text[0] == '{' and text[-1] == '}'
    out=[]; start=1; depth=0
    for i,c in enumerate(text[1:-1],1):
        if c in '{[(': depth+=1
        elif c in '}])': depth-=1
        elif c==',' and depth==0:
            out.append(text[start:i].strip()); start=i+1
        assert depth>=0
    assert depth==0
    out.append(text[start:-1].strip())
    return out

def source_rows(path):
    text=re.sub(r'\(\*.*?\*\)','',path.read_text(encoding='utf-8'),flags=re.S)
    return [split_list(r) for r in split_list(text)]

def parse(text):
    text=text.replace(r'\[Epsilon]', '0').replace('^','**')
    text=' '.join(text.split())
    def walk(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            return FF(node.value)
        if isinstance(node,ast.Name):
            assert node.id in NAMES, node.id
            return NAMES[node.id]
        if isinstance(node,ast.UnaryOp):
            v=walk(node.operand)
            if isinstance(node.op,ast.USub): return -v
            if isinstance(node.op,ast.UAdd): return v
        if isinstance(node,ast.BinOp):
            if isinstance(node.op,ast.Pow):
                e=node.right
                sign=1
                if isinstance(e,ast.UnaryOp) and isinstance(e.op,ast.USub):
                    sign=-1; e=e.operand
                assert isinstance(e,ast.Constant) and type(e.value) is int
                return walk(node.left)**(sign*e.value)
            a,b=walk(node.left),walk(node.right)
            if isinstance(node.op,ast.Add): return a+b
            if isinstance(node.op,ast.Sub): return a-b
            if isinstance(node.op,ast.Mult): return a*b
            if isinstance(node.op,ast.Div): return a/b
        raise ValueError(ast.dump(node))
    return walk(ast.parse(text,mode='eval').body)

def poly_out(p):
    return [[*e,str(c)] for e,c in sorted(p.items())]

def poly_in(p):
    return RR.from_dict({tuple(r[:4]):QQ(r[4]) for r in p})

def enc(x):
    return {'n':poly_out(x.numer),'d':poly_out(x.denom)}

def dec(x):
    return FF.new(poly_in(x['n']),poly_in(x['d']))

def save(name,data):
    OUT.mkdir(exist_ok=True)
    p=OUT/name; tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n',encoding='utf-8')
    tmp.replace(p)

def load_matrix(name):
    data=json.loads((OUT/name).read_text(encoding='utf-8'))
    return [[dec(v) for v in row] for row in data['matrix']]

def degree(p):
    return max((sum(e) for e in p),default=-1)

def mm(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),FF.zero)
             for j in range(len(B[0]))] for i in range(len(A))]

def transpose(A): return list(map(list,zip(*A)))
def diff(A,k): return [[v.diff(XX[k]) for v in row] for row in A]
def sub(A,B): return [[a-b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def add(A,B): return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
