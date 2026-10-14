"""Exact pole-clearance certificate for based rectangular loops; no transport."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
source=folder/'pairing_and_roots.json'
roots=json.loads(source.read_text(encoding='utf-8'))
boxes=[]
for a in roots['rational_finite_poles']:
    boxes.append((str(a),[F(a),F(a)],[F(0),F(0)]))
for i,b in enumerate(roots['q5_boxes'],1):
    boxes.append((f'q5_root_{i}',list(map(F,b['re'])),list(map(F,b['im']))))

def gap(a,b):return max(F(0),a[0]-b[1],b[0]-a[1])
def bounds(p,q):return [[min(p[i],q[i]),max(p[i],q[i])] for i in range(2)]
def clearance(p,q,box):
    assert p[0]==q[0] or p[1]==q[1]
    rect=bounds(p,q)
    return max(gap(rect[0],box[1]),gap(rect[1],box[2]))

base=(F(-64),F(0));r=F(1,4)
records=[]
for a in roots['rational_finite_poles']:
    a=F(a)
    stem=[base,(F(-64),F(1)),(a,F(1)),(a,r)]
    circle=[(a,r),(a-r,r),(a-r,-r),(a+r,-r),(a+r,r),(a,r)]
    twice_area=sum(p[0]*q[1]-q[0]*p[1] for p,q in zip(circle,circle[1:]))
    assert twice_area>0
    interior=([a-r,a+r],[-r,r])
    winding={}
    for name,bx,by in boxes:
        inside=(interior[0][0]<bx[0]<=bx[1]<interior[0][1] and
                interior[1][0]<by[0]<=by[1]<interior[1][1])
        outside=gap(interior[0],bx)>0 or gap(interior[1],by)>0
        assert inside or outside
        winding[name]=1 if inside else 0
    assert [name for name,n in winding.items() if n]==[str(a)]
    vertices=stem+circle[1:]+list(reversed(stem[:-1]))
    certificates=[]
    for p,q in zip(vertices,vertices[1:]):
        distances={name:clearance(p,q,(name,bx,by)) for name,bx,by in boxes}
        d=min(distances.values())
        assert d>=r
        certificates.append({'lower_bound_euclidean_distance':str(d),
                             'nearest_boxes':[k for k,v in distances.items() if v==d]})
    records.append({'center':str(a),'vertices':[[str(v) for v in p] for p in vertices],
                    'segment_clearance':certificates,'winding_numbers':winding})
out={'status':'nine_based_loops_with_exact_clearance_certified',
     'system_sha256':roots['system_sha256'],
     'root_certificate_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'basepoint':['-64','0'],'loops':records,'minimum_path_pole_distance':'1/4',
     'uniform_safe_Taylor_disk_radius':'1/8',
     'all_infinity_distances_unneeded_for_bounded_paths':True,
     'root_boxes_contain_all_finite_poles':True,
     'transport_certified':False,'comparison_to_old_loop_conventions_certified':False,
     'period_initial_conditions_certified':False}
(folder/'paths.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'minimum_clearance':'1/4','loops':len(records)}))
