"""Replay eight independent local calibrations against saved source coefficients."""
from node_frobenius import *


def main():
    p=FOLDER/'marked_lattice/all_node_scale_crosscheck.json';d=json.loads(p.read_text())
    system=FOLDER/'polarized_system.json';assert d['system_sha256']==sha(system)
    data=json.loads(system.read_text());ss=sp.Symbol('ss');D=sp.diag(*map(sp.Rational,D_VALUES))
    K=D.inv()*sp.Matrix(data['inverse_pairing'])*D.inv()
    assert [v['node'] for v in d['nodes']]==POLES[1:]
    for item in d['nodes']:
        a=item['node'];spth=FOLDER/f'node_frobenius/node_{a}.json';seed=json.loads(spth.read_text())
        assert sha(spth)==item['seed_sha256']
        r=list(map(sp.Integer,item['r']));assert [abs(v) for v in r]==[1,2,4,8] and r[3]==8 and sum(r)**2==a
        Q=sp.Matrix(seed['Q'])
        v=[Q*sp.Matrix([sp.Rational(x)+sp.I*sp.Rational(y) for x,y in p]) for p in seed['frobenius_vector_coefficients_exact'][:2]]
        assert v[0][0]==0 and v[1][0]==sp.sympify(item['first_coefficient'])
        K0=K.subs(ss,a);assert (v[0].T*K0*v[0])[0]==0
        norm=sp.simplify(sum(v[1][i]*K0[i,j]*v[0][j]+v[0][i]*K0[i,j]*v[1][j]+sp.I/4*v[0][i]*sp.diff(K[i,j],ss).subs(ss,a)*v[0][j] for i in range(6) for j in range(6)))
        rho2=sp.simplify(4096*sp.I/(v[1][0]**2*a*sp.prod(r)))
        assert norm==sp.sympify(item['norm']) and rho2==sp.sympify(item['rho2_over_pi2'])
        assert sp.simplify(-rho2*norm/2)==sp.Integer(item['nu_over_pi2'])==8
    out={'status':'all_eight_Morse_scale_crosschecks_replayed','certificate_sha256':sha(p),
         'verifier_sha256':sha(__file__),'nodes':POLES[1:],'common_pairing_constant_over_pi_squared':8,
         'floating_point_used':False}
    save_json(p.with_name('all_node_scale_replay.json'),out);print(json.dumps(out))


if __name__=='__main__':main()
