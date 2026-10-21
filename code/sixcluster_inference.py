"""Country inference for a balanced two-phase single-break FE comparison.

Only phase totals, never invented event months, are inputs. Country-period
means identify the coefficient and relevant score subspace for this design.
Not valid for medicine FE models, trends, arbitrary missing country-months,
or covariates unless their reduction to this design has been established.
"""
import itertools
import numpy as np
from scipy.stats import t

def design(d, treated):
    d=np.asarray(d,float); z=np.asarray(treated,float)
    assert len(d)==len(z) and np.isfinite(d).all() and set(z)=={0.,1.}
    X=np.column_stack([np.ones(len(d)),z]);B=np.linalg.inv(X.T@X)
    a=(B@X.T)[1];H=X@B@X.T;M=np.eye(len(d))-H
    b=float(a@d);r=M@d
    return d,z,X,a,H,M,b,r

def weight_support(G,kind):
    support=[-1.,1.] if kind=='rademacher' else [-np.sqrt(1.5),-1.,-np.sqrt(.5),np.sqrt(.5),1.,np.sqrt(1.5)]
    return np.asarray(list(itertools.product(support,repeat=G)))

def wild_pvalue(d,treated,null=0.,kind='webb',pre_months=28,post_months=12,weights=None):
    d,z,X,a,H,M,b,r=design(d,treated);G=len(d);T=pre_months+post_months
    factor=G/(G-1)*(G*T-1)/(G*T-T-G)
    se=np.sqrt(factor*np.sum((a*r)**2));tobs=(b-null)/se
    # Restricted OLS under beta=null, with the common change unrestricted.
    u0=d-null*z;mu=u0.mean();u0=u0-mu
    W=weight_support(G,kind) if weights is None else weights
    ys=null*z+mu+W*u0
    bs=ys@a;rs=ys@M
    ses=np.sqrt(factor*np.sum((rs*a)**2,axis=1))
    numer=bs-null
    fallback=np.zeros_like(numer);nz=np.abs(numer)>=1e-14
    fallback[nz]=np.copysign(np.inf,numer[nz])
    ts=np.divide(numer,ses,out=fallback,where=ses>1e-14)
    # Exhaustive support, so no Monte Carlo +1 correction. Numerical tolerance
    # retains equality on the original and sign-reversed residual draws.
    p=float(np.mean(np.abs(ts)>=abs(tobs)-1e-10))
    return p,float(tobs),len(W)

def wild_acceptance_set(d,treated,kind='webb',pre_months=28,post_months=12,alpha=.05):
    """Numerical inversion; report all detected components, with explicit grid.

    Fine 257-point scan over +/-16 conventional SE, expanding to +/-256 SE
    when necessary; bisection locates each detected crossing to 1e-7 SE.
    Finite-support p functions can be discontinuous; this is a numerical
    acceptance-set diagnostic rather than a claim of analytic endpoints.
    """
    d,z,X,a,H,M,b,r=design(d,treated);G=len(d);T=pre_months+post_months
    factor=G/(G-1)*(G*T-1)/(G*T-T-G)
    se=float(np.sqrt(factor*np.sum((a*r)**2)));W=weight_support(G,kind)
    def p(x):return wild_pvalue(d,z,x,kind,pre_months,post_months,W)[0]
    span=16*se
    for _ in range(5):
        grid=np.linspace(b-span,b+span,257);ps=np.array([p(x) for x in grid]);ok=ps>=alpha
        if not(ok[0] or ok[-1]):break
        span*=2
    components=[];edges=np.flatnonzero(ok[1:]!=ok[:-1]);crossings=[]
    for j in edges:
        lo,hi=grid[j:j+2];state=ok[j]
        for _ in range(32):
            mid=(lo+hi)/2
            if (p(mid)>=alpha)==state:lo=mid
            else:hi=mid
        crossings.append((lo+hi)/2)
    start=-np.inf if ok[0] else None
    for j,x in enumerate(crossings):
        if start is None:start=x
        else:components.append([float(start),float(x)]);start=None
    if start is not None:components.append([float(start),float('inf')])
    return {'components':components,'scan_min':float(grid[0]),'scan_max':float(grid[-1]),'scan_points':257,'alpha':alpha,'endpoint_tolerance':se*2**-27,'warning':'Finite grid may miss acceptance components narrower than grid spacing; endpoints are numerical inversions, not analytic.'}

def infer_changes(changes,treated,pre_months=28,post_months=12,bootstrap=True,invert=False):
    d,z,X,a,H,M,b,r=design(changes,treated);G=len(d);T=pre_months+post_months
    factor=G/(G-1)*(G*T-1)/(G*T-T-G)
    se1=float(np.sqrt(factor*np.sum((a*r)**2)));df1=G-1
    # Genuine CR2 under an identity working covariance, reduced to the
    # country-change score subspace. diag(H)=1/(number in spending arm).
    adj=1/np.sqrt(1-np.diag(H));se2=float(np.sqrt(np.sum((a*adj*r)**2)))
    A=(a*adj)[:,None]*M
    Q=A.T@A;df2=float(np.trace(Q)**2/np.trace(Q@Q))
    out={'estimate':b,'country_clusters':G,'treated_clusters':int(z.sum()),'control_clusters':int(G-z.sum()),'pre_months':pre_months,'post_months':post_months,'CR1_SE':se1,'CR1_df':df1,'CR1_lower95':b-t.ppf(.975,df1)*se1,'CR1_upper95':b+t.ppf(.975,df1)*se1,'CR1_p_two_sided':2*t.sf(abs(b/se1),df1),'CR2_SE':se2,'CR2_Satterthwaite_df':df2,'CR2_lower95':b-t.ppf(.975,df2)*se2,'CR2_upper95':b+t.ppf(.975,df2)*se2,'CR2_p_two_sided':2*t.sf(abs(b/se2),df2),'working_covariance':'Identity country-month covariance; balanced country-change score subspace','status':'valid completed'}
    if bootstrap:
        for kind in ['rademacher','webb']:
            p,tobs,n=wild_pvalue(d,z,kind=kind,pre_months=pre_months,post_months=post_months)
            out.update({f'WCR_{kind}_p_two_sided':p,f'WCR_{kind}_draws':n})
            if invert:out[f'WCR_{kind}_95_acceptance_set']=wild_acceptance_set(d,z,kind,pre_months,post_months)
    return out
