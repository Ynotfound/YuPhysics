"""Explicit horizontal route including all nine upper branch landmarks."""
from refined_transport import *

def main():
    dest=FOLDER/'horizontal_trunk';dest.mkdir(exist_ok=True)
    seedp=FOLDER/'high_trunk_transport.json';seed=json.loads(seedp.read_text())
    rep=json.loads((FOLDER/'high_trunk_replay.json').read_text())
    assert rep['certificate_sha256']==sha(seedp) and seed['end']==['-64','64']
    c=(-64,64);steps=[]
    for a in POLES:
        while c[0]<a:
            end=(min(c[0]+4,a),64);h=ga(end,gn(c));lower,ds=disk_lower(c,h,1)
            steps.append({'start':list(map(str,c)),'end':list(map(str,end)),
                          'step':list(map(str,h)),'disk_radius_upper':str(ab(h)),
                          'pole_distance_lower_on_disk':list(map(str,ds)),
                          'monic_denominator_lower':str(lower),
                          'landmark_pole':a if end[0]==a else None})
            c=end
    out={'status':'horizontal_disks_certified_before_transport','system_sha256':sha(FOLDER/'polarized_system.json'),
         'seed_file':'../high_trunk_transport.json','seed_sha256':sha(seedp),'start':['-64','64'],
         'end':['225','64'],'poles':POLES,'max_step':4,'steps':steps,
         'code_sha256':sha(__file__),'transport_certified':False}
    p=dest/'plan.json'
    if p.exists(): assert json.loads(p.read_text())==out
    else: save_json(p,out)
    print(json.dumps({'status':out['status'],'steps':len(steps),'minimum_clearance':str(min(min(map(F,s['pole_distance_lower_on_disk'])) for s in steps))}))

if __name__=='__main__': main()
