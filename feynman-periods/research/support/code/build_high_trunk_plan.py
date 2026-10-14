"""Save disk-clearance witnesses before any dependent high-trunk integration."""
from gaussian_transport import *

def main():
    points=[(-64,y) for y in [1,9,17,25,33,41,49,57,64]]
    # Includes the future common horizontal trunk; current integration stops
    # at (-64,64). Step radii remain below the exact L-infinity pole distance.
    points += [(x,64) for x in list(range(-56,225,8))+[225]]
    steps=[]
    for c,end in zip(points,points[1:]):
        h=ga(end,gn(c));lower,ds=disk_lower(c,h,1)
        steps.append({'start':list(map(str,c)),'end':list(map(str,end)),
                      'step':list(map(str,h)),'disk_radius_upper':str(ab(h)),
                      'pole_distance_lower_on_disk':list(map(str,ds)),
                      'monic_denominator_lower':str(lower)})
    out={'status':'high_trunk_disks_certified_before_transport','system_sha256':sha(FOLDER/'polarized_system.json'),
         'code_sha256':sha(__file__),'poles':POLES,'path_class_proof':'../../POLARIZED_TRANSPORT.md#三新闭路与原六循环约定的比较',
         'reuse_seed':'first_trunk_step.json','steps':steps,'vertical_steps':8,
         'transports_certified_by_this_file':False}
    (FOLDER/'high_trunk_plan.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'steps':len(steps),'minimum_disk_clearance':min(min(map(F,s['pole_distance_lower_on_disk'])) for s in steps)},default=str))

if __name__=='__main__': main()
