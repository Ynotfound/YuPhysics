"""Two independently factorized four-body phase spaces; numerical estimates, not intervals."""
from pathlib import Path
from itertools import combinations, permutations
import hashlib,json,math,time
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/phi4_four_body'
PARTS=[(0,)+p for p in combinations(range(1,6),2)]

def unit(u,v):
    z=2*u-1; r=np.sqrt(np.maximum(0,1-z*z)); phi=2*np.pi*v
    return np.stack((r*np.cos(phi),r*np.sin(phi),z),axis=-1)

def lamroot(S,a,b):
    return np.sqrt(np.maximum(0,(S-(np.sqrt(a)+np.sqrt(b))**2)*(S-(np.sqrt(a)-np.sqrt(b))**2)))

def boost(k,P,M):
    # Rest to lab, evaluated without cancellation at zero velocity.
    e=k[:,0]; v=k[:,1:]; dot=np.sum(v*P[:,1:],axis=1)
    return np.column_stack(((P[:,0]*e+dot)/M,
           v+P[:,1:]*(e/M+dot/(M*(P[:,0]+M)))[:,None]))

def decay(S,a,b,n):
    M=np.sqrt(S); k=lamroot(S,a,b)/(2*M)
    pa=np.column_stack(((S+a-b)/(2*M),k[:,None]*n))
    pb=np.column_stack(((S+b-a)/(2*M),-k[:,None]*n))
    return pa,pb

def smooth(u):
    return (1-np.cos(np.pi*u))/2, np.pi*np.sin(np.pi*u)/2

def pair_events(s,U):
    n=len(U); one=np.ones(n); ss=s*one
    u,ju=smooth(U[:,0]);v,jv=smooth(U[:,1])
    da=(np.sqrt(s)-2)**2-4; a=4+da*u
    db=(np.sqrt(s)-np.sqrt(a))**2-4;b=4+db*v
    Q,R=decay(ss,a,b,unit(U[:,2],U[:,3]))
    k1,k2=decay(a,one,one,unit(U[:,4],U[:,5]))
    k3,k4=decay(b,one,one,unit(U[:,6],U[:,7]))
    ks=np.stack((boost(k1,Q,np.sqrt(a)),boost(k2,Q,np.sqrt(a)),
                 boost(k3,R,np.sqrt(b)),boost(k4,R,np.sqrt(b))),axis=1)
    weight=da*db*ju*jv*lamroot(ss,a,b)/s*np.sqrt(1-4/a)*np.sqrt(1-4/b)/(2048*np.pi**5)
    return ks,weight

def sequential_events(s,U):
    n=len(U);one=np.ones(n);ss=s*one
    u,ju=smooth(U[:,0]);v,jv=smooth(U[:,1])
    da=(np.sqrt(s)-1)**2-9;a=9+da*u
    db=(np.sqrt(a)-1)**2-4;b=4+db*v
    Q,k4=decay(ss,a,one,unit(U[:,2],U[:,3]))
    R,k3=decay(a,b,one,unit(U[:,4],U[:,5]))
    k1,k2=decay(b,one,one,unit(U[:,6],U[:,7]))
    k1=boost(k1,R,np.sqrt(b));k2=boost(k2,R,np.sqrt(b))
    ks=np.stack((boost(k1,Q,np.sqrt(a)),boost(k2,Q,np.sqrt(a)),
                 boost(k3,Q,np.sqrt(a)),k4),axis=1)
    weight=da*db*ju*jv*lamroot(ss,a,one)/s*lamroot(a,b,one)/a*np.sqrt(1-4/b)/(2048*np.pi**5)
    return ks,weight

def amplitude(s,ks):
    n=len(ks);p=np.sqrt(s/4-1);E=np.sqrt(s)/2
    incoming=np.broadcast_to(np.array([[-E,0,0,-p],[-E,0,0,p]]),(n,2,4))
    qs=np.concatenate((incoming,ks),axis=1)
    den=np.stack([np.sum(qs[:,A,:],axis=1) for A in PARTS],axis=1)
    den=den[:,:,0]**2-np.sum(den[:,:,1:]**2,axis=2)-1
    return np.sum(1/den,axis=1),den

def invariants(s,ks):
    total=np.sum(ks,axis=1);total[:,0]-=np.sqrt(s)
    shell=ks[:,:,0]**2-np.sum(ks[:,:,1:]**2,axis=2)-1
    return float(np.max(np.abs(total))),float(np.max(np.abs(shell)))

def run(s,factory,powers=(14,17,19),reps=8):
    rows=[]
    for power in powers:
        stats=[];err=[0.,0.];mins=np.inf;maxs=-np.inf
        for j in range(reps):
            sampler=qmc.Sobol(d=8,scramble=True,seed=90107+7919*j+(factory==sequential_events)*107)
            totals=np.zeros(3)
            # Batch by powers of two; no points discarded.
            blockpower=min(power,15)
            for _ in range(2**(power-blockpower)):
                U=sampler.random(2**blockpower)
                ks,w=factory(s,U);S,den=amplitude(s,ks)
                totals += [np.sum(w),np.sum(w*S*S),np.sum(w*np.sum(1/den**2,axis=1))]
                a,b=invariants(s,ks);err=[max(a,err[0]),max(b,err[1])]
                mins=min(mins,float(np.min(den)));maxs=max(maxs,float(np.max(den)))
            stats.append(totals/2**power)
        arr=np.array(stats);mean=arr.mean(0);se=arr.std(axis=0,ddof=1)/np.sqrt(reps)
        fluxsym=48*s*np.sqrt(1-4/s)
        rows.append(dict(power=power,replicates=reps,mean_phase_volume=mean[0],
                         sigma_over_g4=mean[1]/fluxsym,replicate_standard_error=se[1]/fluxsym,
                         phase_replicate_standard_error=se[0],
                         diagonal_only_sigma=mean[2]/fluxsym,
                         replicate_sigma=(arr[:,1]/fluxsym).tolist(),
                         conservation_max=err[0],shell_max=err[1],denominator_range=[mins,maxs]))
        print(s,factory.__name__,power,rows[-1]['sigma_over_g4'],rows[-1]['replicate_standard_error'],flush=True)
    return rows

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data=dict(status='numerical_estimates_not_certified_intervals',mass=1,identical_factorial=24,
              partitions=PARTS,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              rows=[])
    for s in [20.,32.]:
        for name,fun in [('pair',pair_events),('sequential',sequential_events)]:
            row=dict(s_over_m2=s,factorization=name,levels=run(s,fun))
            data['rows'].append(row)
            (OUT/'production.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    print('done',flush=True)
if __name__=='__main__':main()

