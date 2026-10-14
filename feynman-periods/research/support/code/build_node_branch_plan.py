"""Save safe ordinary disks connecting certified node seeds to the upper trunk."""
from refined_transport import *


def disk(c,h):
    lower,distances=disk_lower(c,h,1)
    return {'center':list(map(str,c)),'step':list(map(str,h)),
            'radius':str(ab(h)),'pole_distance_lower':list(map(str,distances)),
            'monic_denominator_lower':str(lower)}


def main():
    local=FOLDER/'node_frobenius';horizontal=FOLDER/'horizontal_trunk'
    replay=json.loads((local/'replay.json').read_text())
    landmarks=json.loads((horizontal/'landmarks.json').read_text())
    marks={m['pole']:m for m in landmarks['landmarks']}
    branches=[]
    for entry in replay['certificates']:
        a=entry['node'];p=local/entry['file'];seed=json.loads(p.read_text())
        assert entry['sha256']==sha(p)
        assert seed['endpoint']==[str(a),'1/4'] and marks[a]['point']==[str(a),'64']
        # Eight equal steps on each of eight dyadic height bands.
        y=F(64);steps=[]
        for band in range(8):
            dy=y/16
            for j in range(8):
                c=(F(a),y);h=(F(0),-dy);end=ga(c,h)
                steps.append({'start':list(map(str,c)),'end':list(map(str,end)),
                              'down':disk(c,h),'up':disk(end,gn(h))})
                y-=dy
        assert y==F(1,4) and len(steps)==64
        branches.append({'node':a,'upper':[str(a),'64'],'local':[str(a),'1/4'],
                         'local_certificate':entry['file'],'local_sha256':sha(p),
                         'steps_top_to_bottom':steps})
    out={'status':'eight_bidirectional_node_branch_disk_plans_certified',
         'system_sha256':sha(FOLDER/'polarized_system.json'),
         'landmarks_sha256':sha(horizontal/'landmarks.json'),
         'local_replay_sha256':sha(local/'replay.json'),'code_sha256':sha(__file__),
         'poles':POLES,'branches':branches,'ordinary_segment_count':512,
         'directed_disk_count':1024,'descent_minimum_clearance':'1/4',
         'ascent_minimum_clearance':'7/32','ordinary_transport_certified':False,
         'zero_branch_included':False,
         'intended_use':'Transport the certified local vanishing line upwards; then apply the certified horizontal inverse.'}
    target=local/'branch_plan.json'
    if target.exists(): assert json.loads(target.read_text())==out
    else: save_json(target,out)
    print(json.dumps({k:v for k,v in out.items() if k!='branches'}))


if __name__=='__main__': main()
