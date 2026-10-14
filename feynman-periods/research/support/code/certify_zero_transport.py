"""Transport a certified zero Frobenius frame and recover the based monodromy."""
from zero_frobenius import *
from certify_high_trunk import proposal
from refined_transport import compose_fine
from node_branch_transport import propagate_vector,reflection_at_vector
import time


def matrix_switch(C,e,T,bits):
    S=1<<bits;exact=mm(T,[gs(F(1,S),v) for v in C])
    new=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact,new)])
    return new,radius_up(mnorm(T)*e+rd,bits),rd


def eval_rational_matrix(M,point):
    s=sp.Symbol('ss');z=sp.Rational(point[0])+sp.I*sp.Rational(point[1])
    A=M.subs(s,z).applyfunc(sp.expand)
    return [(rat(sp.re(v)),rat(sp.im(v))) for v in A]


def pi_certificate(terms=100):
    def arctan_bounds(q):
        x=F(1,q);partial=sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms))
        next_term=(-1)**terms*x**(2*terms+1)/F(2*terms+1)
        return min(partial,partial+next_term),max(partial,partial+next_term)
    l5,u5=arctan_bounds(5);l239,u239=arctan_bounds(239)
    low=16*l5-4*u239;high=16*u5-4*l239
    return {'formula':'pi=16*atan(1/5)-4*atan(1/239)','terms':terms,
            'lower':str(low),'upper':str(high),'center':str((low+high)/2),'radius':str((high-low)/2)}


def build_plan(dest,seedp):
    y=F(1,4);steps=[]
    for band in range(8):
        dy=y/8
        for j in range(8):
            c=(F(0),y);h=(F(0),dy);lower,ds=disk_lower(c,h,1)
            steps.append({'start':list(map(str,c)),'end':['0',str(y+dy)],'step':['0',str(dy)],
                          'pole_distance_lower':list(map(str,ds)),'monic_denominator_lower':str(lower)})
            y+=dy
    assert y==64 and len(steps)==64
    out={'status':'zero_ascent_disks_certified','system_sha256':sha(FOLDER/'polarized_system.json'),
         'seed_sha256':sha(seedp),'code_sha256':sha(__file__),'start':['0','1/4'],'end':['0','64'],
         'switch_to_D_height':'2','steps':steps,'transport_certified':False}
    p=dest/'plan.json'
    if p.exists(): assert json.loads(p.read_text())==out
    else: save_json(p,out)
    return out,p


def main():
    dest=FOLDER/'zero_monodromy';seedp=dest/'frobenius.json';seed=json.loads(seedp.read_text())
    assert json.loads((dest/'frobenius_replay.json').read_text())['certificate_sha256']==sha(seedp)
    plan,planp=build_plan(dest,seedp);bits=seed['bits'];S=1<<bits
    zero_source=(list(map(F,seed['full_source_q'])),[list(map(F,M)) for M in seed['full_source_N']])
    source=source_polynomials();toD=sp.Matrix(seed['to_D'])
    C=decode(seed['endpoint_center_gaussian_integer']);e=F(seed['endpoint_error_upper']);frame='zero_Q'
    summary={'status':'partial_zero_Frobenius_transport','bits':bits,
             'system_sha256':sha(FOLDER/'polarized_system.json'),'seed_sha256':sha(seedp),
             'plan_sha256':sha(planp),'code_sha256':sha(__file__),
             'start':seed['endpoint'],'segments':[],'absolute_periods_certified':False}
    started=time.time()
    for ix,step in enumerate(plan['steps']):
        tick=time.time();c=tuple(map(F,step['start']));h=tuple(map(F,step['step']))
        switched=False;switch_error=F(0);T=None
        if frame=='zero_Q' and c[1]>=2:
            T=eval_rational_matrix(toD,c);C,e,switch_error=matrix_switch(C,e,T,bits)
            frame='D';switched=True
        src=zero_source if frame=='zero_Q' else source
        q,N,L=local_polynomials(src,c,h);lower,ds=disk_lower(c,h,L)
        assert list(map(str,ds))==step['pole_distance_lower']
        p=dest/f'step_{ix+1:03d}.json'
        if p.exists():
            old=json.loads(p.read_text());assert old['frame']==frame and old['plan_sha256']==sha(planp)
            P=[decode(M) for M in old['polynomial_coefficients_gaussian_integer']]
        else: P=proposal(q,N,96,bits)
        M,bd=bounds(q,N,P,bits,lower);le=radius_up(bd['error_norm_upper'],bits)
        assert le<F(2)**-170, f'zero step {ix+1} requires refinement'
        new,ne,rd=compose_fine(M,le,C,e,bits)
        entry={'status':'zero_frame_step_certified','index':ix+1,'bits':bits,'order':len(P)-1,
               'plan_sha256':sha(planp),'system_sha256':summary['system_sha256'],'frame':frame,
               'start':step['start'],'end':step['end'],'switched_before_step':switched,
               'switch_matrix_exact':None if T is None else encode(T),'switch_rounding_error':str(switch_error),
               'input_center_gaussian_integer':encode(C),'input_error_upper':str(e),
               'common_q_N_scale':str(L),'q_coefficients_gaussian_integer':encode(q),
               'N_coefficients_gaussian_integer':[encode(A) for A in N],
               'polynomial_coefficients_gaussian_integer':[encode(A) for A in P],
               'bounds':{k:str(v) for k,v in bd.items()},'local_matrix_error_upper':str(le),
               'transition_center_gaussian_integer':encode(M),
               'output_center_gaussian_integer':encode(new),'output_error_upper':str(ne),
               'matrix_rounding_error':str(rd)}
        if p.exists(): assert old==entry
        else: save_json(p,entry)
        C,e=new,ne
        summary['segments'].append({'file':p.name,'sha256':sha(p),'end':step['end']})
        summary.update(end=step['end'],frame=frame,center_gaussian_integer=encode(C),error_upper=str(e))
        save_json(dest/'transport.json',summary)
        if ix==0 or (ix+1)%8==0 or switched:
            print(json.dumps({'zero_step':ix+1,'height':step['end'][1],'frame':frame,
                              'M_display':float(bd['matrix_norm_upper']),
                              'error_display':float(e),'seconds':round(time.time()-tick,2)}),flush=True)
    hp=FOLDER/'horizontal_trunk/landmarks.json';hor=json.loads(hp.read_text())
    item=next(x for x in hor['landmarks'] if x['pole']==0)
    factor=1<<(bits-hor['bits']);H=[gs(factor,v) for v in decode(item['high_inverse']['center_gaussian_integer'])]
    he=F(item['high_inverse']['error_upper'])
    W,we,rd=compose_fine(H,he,C,e,bits);winv=inverse(W,we,bits)
    Vi=decode(winv['center_gaussian_integer']);ve=F(winv['error_upper'])
    Nint=[gs(S,v) for v in real_matrix(sp.Matrix(seed['nilpotent']))]
    WN,wne,rd1=compose_fine(W,we,Nint,F(0),bits)
    A,ae,rd2=compose_fine(WN,wne,Vi,ve,bits)
    picert=pi_certificate();pc=F(picert['center']);pe=F(picert['radius'])
    Ac=[gs(F(1,S),v) for v in A]
    Texact=ma(mi(),ms((F(0),2*pc),Ac))
    T=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in Texact]
    trd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(Texact,T)])
    te=radius_up(2*((abs(pc)+pe)*ae+pe*mnorm(Ac))+trd,bits)
    n1p=FOLDER/'node_branches/node_1/transport.json';n1=json.loads(n1p.read_text())
    v1=[gs(1<<(bits-n1['bits']),v) for v in decode(n1['based_center_gaussian_integer'])]
    e1=F(n1['based_error_upper'])
    v6,e6,rd6=propagate_vector(T,te,v1,e1,bits)
    K=pairing_at(-64,(F(0),F(64)))
    ref6=reflection_at_vector([gs(F(1,S),v) for v in v6],e6,K,bits)
    summary.update(status='zero_based_monodromy_and_sixth_reflection_certified',
        landmarks_sha256=sha(hp),node_one_sha256=sha(n1p),based_point=['-64','64'],
        based_frame_center_gaussian_integer=encode(W),based_frame_error_upper=str(we),
        rebasing_rounding_error=str(rd),based_frame_inverse=winv,
        conjugated_nilpotent_center_gaussian_integer=encode(A),conjugated_nilpotent_error_upper=str(ae),
        nilpotent_first_rounding_error=str(rd1),nilpotent_second_rounding_error=str(rd2),
        pi_certificate=picert,monodromy_center_gaussian_integer=encode(T),monodromy_error_upper=str(te),
        monodromy_rounding_error=str(trd),sixth_vector_center_gaussian_integer=encode(v6),
        sixth_vector_error_upper=str(e6),sixth_vector_rounding_error=str(rd6),
        based_pairing_exact=encode(K),sixth_reflection=ref6)
    save_json(dest/'transport.json',summary)
    print(json.dumps({'status':summary['status'],'T0_error_dyadic':str(dyadic_upper(te)),
                      'sixth_reflection_error_dyadic':str(dyadic_upper(F(ref6['reflection_error_upper']))),
                      'seconds':round(time.time()-started,2)}),flush=True)


if __name__=='__main__': main()
