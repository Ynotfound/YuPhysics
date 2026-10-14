"""Fine outward radius rounding without changing the archived first backend."""
from gaussian_transport import *

def radius_up(x,bits):
    assert x>=0
    S=1<<bits
    return F(ceil(F(x)*S),S)

def compose_fine(L,le,C,ce,bits):
    S=1<<bits;product=mm(L,C)
    new=[(dyadic_round(F(a,S*S),bits),dyadic_round(F(b,S*S),bits)) for a,b in product]
    rounding=max(sum(F(ab(ga(product[6*i+j],gn(gs(S,new[6*i+j])))),S*S)
                     for j in range(6)) for i in range(6))
    raw=(F(mnorm(L),S)+le)*ce+le*F(mnorm(C),S)+rounding
    upper=radius_up(raw,bits)
    assert raw<=upper<raw+F(1,S)
    return new,upper,rounding

def save_json(path,value):
    # Same-directory atomic replacement only for this driver's own summaries.
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)
