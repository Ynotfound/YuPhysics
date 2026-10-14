"""Independent rational check of both directions of each saved node branch."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/banana_pure_ray'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def check(d):
    poles=[0,1,9,25,49,81,121,169,225];assert d['poles']==poles
    assert [b['node'] for b in d['branches']]==poles[1:]
    minima={'down':None,'up':None};count=0
    for branch in d['branches']:
        a=branch['node'];current=[F(a),F(64)]
        assert branch['upper']==[str(a),'64'] and branch['local']==[str(a),'1/4']
        local=DEST/'node_frobenius'/branch['local_certificate']
        assert sha(local)==branch['local_sha256']
        assert json.loads(local.read_text())['endpoint']==branch['local']
        for step in branch['steps_top_to_bottom']:
            start=list(map(F,step['start']));end=list(map(F,step['end']))
            assert current==start and start[0]==end[0]==a and 0<end[1]<start[1]
            for direction,c,e in [('down',start,end),('up',end,start)]:
                witness=step[direction];h=[e[i]-c[i] for i in range(2)]
                assert list(map(F,witness['center']))==c
                assert list(map(F,witness['step']))==h
                radius=abs(h[0])+abs(h[1]);assert radius==F(witness['radius'])
                lower=F(1);distances=[]
                for p in poles:
                    distance=max(abs(c[0]-p),abs(c[1]))-radius
                    assert distance>0
                    # Reverse triangle inequality: the Euclidean distance on
                    # the entire closed disk is at least this rational value.
                    assert max(abs(c[0]-p),abs(c[1]))**2 <= (c[0]-p)**2+c[1]**2
                    lower*=distance;distances.append(distance)
                assert distances==list(map(F,witness['pole_distance_lower']))
                assert lower==F(witness['monic_denominator_lower'])
                minima[direction]=min([min(distances)]+([] if minima[direction] is None else [minima[direction]]))
            current=end;count+=1
        assert current==[F(a),F(1,4)]
    assert count==d['ordinary_segment_count']==512 and d['directed_disk_count']==2*count
    assert minima=={'down':F(d['descent_minimum_clearance']),'up':F(d['ascent_minimum_clearance'])}
    return minima


def main():
    path=DEST/'node_frobenius/branch_plan.json';d=json.loads(path.read_text())
    assert d['system_sha256']==sha(DEST/'polarized_system.json')
    assert d['landmarks_sha256']==sha(DEST/'horizontal_trunk/landmarks.json')
    assert d['local_replay_sha256']==sha(DEST/'node_frobenius/replay.json')
    assert d['code_sha256']==sha(ROOT/'code/build_node_branch_plan.py')
    minima=check(d)
    bad=json.loads(json.dumps(d));bad['branches'][0]['steps_top_to_bottom'][-1]['end'][1]='0'
    rejected=False
    try: check(bad)
    except AssertionError: rejected=True
    assert rejected
    out={'status':'eight_bidirectional_branch_plans_independently_replayed',
         'plan_sha256':sha(path),'verifier_sha256':sha(Path(__file__)),
         'ordinary_segments':512,'directed_disks':1024,
         'minimum_clearance':{k:str(v) for k,v in minima.items()},
         'tampered_pole_endpoint_rejected':rejected,
         'ordinary_transport_certified':False,'common_basepoint_reflections_certified':False}
    target=path.with_name('branch_plan_replay.json')
    target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))


if __name__=='__main__': main()
