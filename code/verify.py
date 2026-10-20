#!/usr/bin/env python3
"""Exact rational replay: no third-party dependencies and no network access.

Checks the explicit lower-bound family by actually applying its two proximal
maps and gradient, and checks the general slack identity algebra independently
on rational arrays. The latter checks algebra, not convexity interpolation.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import random
import sys


class VerificationError(RuntimeError):
    """Raised when a mathematical replay condition fails."""


def require(condition, message):
    """Fail closed even when Python is run with optimization enabled."""
    if not condition:
        raise VerificationError(message)

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def scale(t,a): return tuple(t*x for x in a)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),Q(0))
def norm2(a): return dot(a,a)
def mean(a): return scale(Q(1,len(a)),tuple(map(sum,zip(*a))))
def soft(x,t): return max(x-t,Q(0)) if x>=0 else min(x+t,Q(0))
def clip(x,t): return min(max(x,-t),t)

def replay_family(K,M,alpha,theta):
    M,alpha,theta=Q(M),Q(alpha),Q(theta)
    horizon=Q(K+1)-theta
    weight=Q(K+2)
    P=(Q(1),M)
    R=weight+M+1
    z=(Q(0),Q(0))
    x0=(horizon/2,Q(0)); u0=scale(1/alpha,x0)
    x,u=x0,u0
    xs=[]; us=[]
    def unit_grad(w):
        t=theta*(M*w[0]-w[1])/(M*M+1)
        return (M*t,-t)
    for k in range(1,K+1):
        # prox_{g*/alpha}(u+x/alpha): shifted box projection
        w=sub(add(u,scale(1/alpha,x)),scale(1/alpha,P))
        unew=tuple(clip(t,R/alpha) for t in w)
        p=add(x,scale(alpha,sub(u,unew)))
        require(p == P, f"unexpected p at K={K}, M={M}, alpha={alpha}, k={k}: {p}")
        require(unit_grad(p) == z, f"gradient at P is nonzero: {unit_grad(p)}")
        w=sub(sub(p,scale(alpha,unew)),unit_grad(p))
        xnew=(soft(w[0],weight),soft(w[1],M))
        require(xnew == (Q(0), M), f"unexpected x at {(K, M, alpha, k)}: {xnew}")
        require(unew == ((horizon-k)/alpha, -M/alpha), f"unexpected u at step {k}: {unew}")
        require(all(abs(t) < R/alpha for t in unew), f"dual iterate hit box boundary at step {k}")
        x,u=xnew,unew
        xs.append(x);us.append(u)
    xb,ub=mean(xs),mean(us)
    f=(weight*abs(xb[0])+M*abs(xb[1]))/alpha
    h=theta*(M*xb[0]-xb[1])**2/(2*alpha*(M*M+1))
    gstar=dot(P,ub)
    gap=f+h+gstar
    D=norm2(x0)+alpha*alpha*norm2(u0)
    ratio=alpha*gap/D
    expected=1/horizon-theta/(horizon*horizon*(M*M+1))
    require(ratio == expected, f"ratio mismatch: {ratio} != {expected}")
    require(D == horizon*horizon/2, f"D mismatch: {D}")
    require(gap == (horizon/2-theta/(2*(M*M+1)))/alpha, f"gap mismatch: {gap}")
    # Hessian is theta times an orthogonal rank-one projection.
    H=((theta*M*M/(M*M+1),-theta*M/(M*M+1)),(-theta*M/(M*M+1),theta/(M*M+1)))
    H2=tuple(tuple(sum(H[i][r]*H[r][j] for r in range(2)) for j in range(2)) for i in range(2))
    require(H2 == tuple(tuple(theta*t for t in row) for row in H), "Hessian projector identity failed")
    require(H[0][0]+H[1][1] == theta, "Hessian trace mismatch")
    require(H[0][0]*H[1][1]-H[0][1]*H[1][0] == 0, "Hessian determinant mismatch")
    ustar=(-weight/alpha,-M/alpha)
    require(all(abs(t) < R/alpha for t in ustar), "saddle dual point is not interior")
    require(add((weight/alpha, M/alpha), ustar) == z, "saddle balance failed")
    return dict(K=K,M=str(M),alpha=str(alpha),theta=str(theta),D=str(D),gap=str(gap),ratio=str(ratio),shortfall=str(1/horizon-ratio))

def check_identity(K,seed,theta,d=3):
    theta=Q(theta)
    horizon=Q(K+1)-theta
    rng=random.Random(seed)
    def scalar():return Q(rng.randrange(-9,10),rng.randrange(1,6))
    def vec():return tuple(scalar() for _ in range(d))
    x0,u0,x,u=vec(),vec(),vec(),vec()
    xs=[vec() for _ in range(K)]
    us=[vec() for _ in range(K)]
    bs=[vec() for _ in range(K)]
    ps=[sub(add(xx,uu),un) for xx,uu,un in zip([x0]+xs[:-1],[u0]+us[:-1],us)]
    aa=[sub(sub(sub(p,uk),bk),xk) for p,uk,bk,xk in zip(ps,us,bs,xs)]
    As=[add(ak,uk) for ak,uk in zip(aa,us)]
    vs=[add(A,b) for A,b in zip(As,bs)]
    xb,ub,Ab,bb=mean(xs),mean(us),mean(As),mean(bs)
    S=scale(Q(K),mean(vs))
    hb,hx=vec(),vec()
    fb,fx,gb,gu,hbv,hxv=(scalar() for _ in range(6))
    hp=[scalar() for _ in range(K)]
    gap=fb+hbv+dot(u,xb)-gu-fx-hxv-dot(ub,x)+gb
    D=norm2(sub(x0,x))+norm2(sub(u0,u))
    Jf=fb-fx-sum(dot(a,sub(xk,x)) for a,xk in zip(aa,xs))/K
    Jg=gb-gu-sum(dot(p,sub(uk,u)) for p,uk in zip(ps,us))/K
    Jh=sum(hbv-hpk+dot(hb,sub(p,xb))+norm2(sub(hb,b))/(2*theta)
           +hpk-hxv+dot(b,sub(x,p))+norm2(sub(hx,b))/(2*theta)
           for p,b,hpk in zip(ps,bs,hp))/K
    vbar=mean(vs)
    Sh=(norm2(sub(add(hb,scale(theta,vbar)),bb))+norm2(sub(hx,bb)))/(2*theta)
    centre=scale(horizon/2,vbar)
    S1=(norm2(sub(sub(x0,x),centre))+norm2(sub(sub(u0,u),centre)))/horizon
    S2=sum(norm2(sub(A,Ab))/2+(1/theta-Q(1,2))*norm2(sub(b,bb)) for A,b in zip(As,bs))/K
    residual=gap-D/horizon-(Jf+Jg+Jh-Sh-S1-S2)
    require(residual == 0, f"identity residual at K={K}, seed={seed}: {residual}")
    return str(residual)

def check_one_step(theta):
    theta=Q(theta)
    horizon=2-theta
    d=theta-1
    shift=max(d,Q(0))
    def f(t):return max(Q(0),d*t+shift)
    x0=u0=horizon/2
    u1=clip(u0+x0,Q(4))
    p1=x0+u0-u1
    prox_input=p1-u1-theta*p1
    candidates=[prox_input,prox_input-d]
    if d:candidates.append(-shift/d)
    values=[f(t)+(t-prox_input)**2/2 for t in candidates]
    xmin=candidates[values.index(min(values))]
    require(u1 == horizon and p1 == 0 and xmin == -1, "one-step trajectory mismatch")
    gap=f(xmin)-f(0)+theta*xmin*xmin/2
    D=x0*x0+u0*u0
    require(gap == horizon/2 and gap/D == 1/horizon, "one-step equality mismatch")
    return dict(theta=str(theta),gap=str(gap),D=str(D),ratio=str(gap/D))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--negative-control",
        choices=("trajectory", "constant"),
        help="deliberately test a false trajectory or constant (must exit nonzero)",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.negative_control == "trajectory":
        case = replay_family(2, Q(1), Q(1), Q(1))
        require(case["gap"] == "999", "deliberately false trajectory-gap control")
    if args.negative_control == "constant":
        require(Q(1, 2) == Q(1, 3), "deliberately false sharp-constant control")

    thetas=[Q(0),Q(1,4),Q(1),Q(3,2),Q(19,10),Q(199,100)]
    cases=[replay_family(K,M,a,theta) for K in range(1,11)
           for M in [Q(1,3),Q(1),Q(10),Q(100)]
           for a in [Q(1,3),Q(1),Q(7,2)] for theta in thetas]
    identities=[check_identity(K,1000*K+j,theta) for K in range(1,13)
                for j in range(5) for theta in thetas[1:]]
    one_steps=[check_one_step(theta) for theta in thetas]
    endpoint_cases=[replay_family(K,M,Q(1),Q(2)) for K in range(2,11)
                    for M in [Q(1,3),Q(1),Q(10),Q(100)]]
    endpoint_identities=[check_identity(K,2000*K+j,Q(2)) for K in range(2,11)
                         for j in range(5)]
    report={
        'arithmetic':'Python fractions.Fraction; exact rational arithmetic',
        'validation':'explicit VerificationError checks; active under python -O',
        'python_optimization_flag':sys.flags.optimize,
        'family_cases':len(cases),
        'family_splitting_iterations_checked':sum(c['K'] for c in cases),
        'family_individual_proximal_evaluations':2*sum(c['K'] for c in cases),
        'family_formula':'alpha*gap/D = 1/a - theta/(a^2*(1+M^2)), a=K+1-theta',
        'identity_cases':len(identities),
        'identity_base_array_samples':12*5,
        'identity_positive_theta_values_per_sample':len(thetas)-1,
        'identity_sampling_note':'60 seeded rational arrays; each checked at 5 positive theta values, giving 300 identity cases',
        'identity_all_residuals_zero':all(x=='0' for x in identities),
        'all_explicit_checks_passed':True,
        'one_step_attainment_cases':one_steps,
        'endpoint_theta_2':{
            'scope':'K>=2; sharp upper bounds and suprema only; no attainment classification',
            'family_cases':len(endpoint_cases),
            'identity_cases':len(endpoint_identities),
            'all_exact_checks_passed':True,
        },
        'cases':cases,
    }
    out=args.output or (Path(__file__).resolve().parents[1]/'results'/'davis_yin_exact.json')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='cases'},indent=2))

if __name__=='__main__':main()
