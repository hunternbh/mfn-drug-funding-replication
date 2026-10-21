"""Controlled projections from final country-month tables; no scraped inputs."""
import numpy as np,pandas as pd
from scipy.stats import t

def X(pos,method):
    pos=np.asarray(pos,float);z=[np.ones(len(pos))]
    if 'linear' in method:z.append(pos/12)
    if 'harmonic' in method:z.extend([np.sin(2*np.pi*pos/12),np.cos(2*np.pi*pos/12)])
    return np.column_stack(z)
def project(y,trainpos,targetpos,method,B,seed):
    xx=X(trainpos,method);zz=X(targetpos,method);pinv=np.linalg.pinv(xx)
    beta=pinv@y;resid=y-xx@beta;resid-=resid.mean();n=len(y);rng=np.random.default_rng(seed)
    def draw(length):
        starts=rng.integers(0,n,size=(B,int(np.ceil(length/12))))
        ix=(starts[:,:,None]+np.arange(12)[None,None,:])%n
        return resid[ix.reshape(B,-1)[:,:length]]
    beta_draw=(xx@beta+draw(n))@pinv.T
    means=(beta_draw@zz.T+draw(len(targetpos))).mean(axis=1)
    return zz@beta,means
def run(root,output):
    d=pd.read_csv(root/'data/long_history_country_month.csv',parse_dates=['month'])
    cc=['Brazil','Australia','Switzerland'];methods=['constant_gap','linear_gap_trend','annual_harmonic','linear_plus_harmonic']
    final=[];validation=[];did=[];annual=[];leave=[]
    for oi,outcome in enumerate(['aggregate','exposure']):
        mat=d.pivot(index='country',columns='month',values=outcome).reindex(cc).to_numpy()
        assert mat.shape==(3,73) and np.isfinite(mat).all()
        gap=mat[1:].mean(axis=0)-mat[0];pre=gap[:60];actual=gap[61:].mean()
        did.append(dict(outcome=outcome,DiD=actual-pre.mean()))
        for j,n in enumerate([24,36,47]):
            for mi,method in enumerate(methods):
                target=np.arange(n+1,n+13);f,z=project(gap[:n],np.arange(n),target,method,2000,1000+oi*100+mi*10+j)
                lo,hi=np.quantile(z,[.025,.975]);a=gap[target].mean();error=a-f.mean()
                validation.append(dict(outcome=outcome,method=method,pre_months=n,actual=a,expected=f.mean(),error=error,absolute_error=abs(error),predictive_lower95=lo,predictive_upper95=hi,covered95=bool(lo<=a<=hi)))
        v=pd.DataFrame(validation);v=v[v.outcome.eq(outcome)];scores=v.groupby('method').absolute_error.mean().reindex(methods);selected=scores.idxmin()
        for mi,method in enumerate(methods):
            f,z=project(pre,np.arange(60),np.arange(61,73),method,5000,20261021+oi*100+mi)
            lo,hi=np.quantile(z,[.025,.975])
            final.append(dict(outcome=outcome,method=method,selected_by_preperiod_validation=(method==selected),predicted_gap=f.mean(),observed_gap=actual,average_projection_error=actual-f.mean(),predictive_gap_lower95=lo,predictive_gap_upper95=hi,predictive_effect_lower95=actual-hi,predictive_effect_upper95=actual-lo,validation_MAE=scores.loc[method]))
        blocks=pre.reshape(5,12).mean(axis=1);half=t.ppf(.975,4)*blocks.std(ddof=1)*np.sqrt(1+1/5)
        annual.append(dict(outcome=outcome,effect=actual-pre.mean(),predictive_effect_lower95=actual-pre.mean()-half,predictive_effect_upper95=actual-pre.mean()+half))
        for j in range(5):leave.append(dict(outcome=outcome,omitted_block=f'May{2020+j}-April{2021+j}',estimate=actual-np.delete(blocks,j).mean()))
    tables={'all_controlled_projection_results.csv':pd.DataFrame(final),'prepolicy_projection_validation.csv':pd.DataFrame(validation),'five_year_observed_DiD.csv':pd.DataFrame(did),'annual_block_prediction_sensitivity.csv':pd.DataFrame(annual),'leave_one_preyear_out.csv':pd.DataFrame(leave)}
    checks=[]
    for name,actual in tables.items():
        expected=pd.read_csv(root/'expected'/name)
        keys=[x for x in ['outcome','method','pre_months','omitted_block'] if x in actual and x in expected]
        z=actual.merge(expected,on=keys,suffixes=('_new','_reference'),validate='one_to_one')
        assert len(z)==len(expected)
        for c in actual.columns:
            if c in keys or not pd.api.types.is_numeric_dtype(actual[c]) or pd.api.types.is_bool_dtype(actual[c]):continue
            assert np.allclose(z[c+'_new'],z[c+'_reference'],atol=1e-8,rtol=1e-7),name+':'+c
        actual.to_csv(output/name,index=False);checks.append(dict(table=name,rows=len(actual)))
    return checks
