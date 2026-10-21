"""Offline replication from final aggregate tables, not publisher archives.
Run: python reproduce.py. Results are written only to build/.
"""
from pathlib import Path
import sys,json,hashlib,platform,subprocess
import numpy as np,pandas as pd
from scipy.stats import t

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'build';TABLES=OUT/'tables'
sys.path.insert(0,str(ROOT/'code'))
from sixcluster_inference import infer_changes
from projections import run as projections

def main():
    TABLES.mkdir(parents=True,exist_ok=True)
    for path,digest in json.loads((ROOT/'data/checksums.json').read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    print('Final-table checksums verified.',flush=True)
    definitions=pd.read_csv(ROOT/'data/alternative_specification_definitions.csv')
    country=pd.read_csv(ROOT/'data/alternative_country_period_outcomes.csv')
    rows=[]
    for spec,q in country.groupby('spec_id',sort=False):
        m=definitions[definitions.spec_id.eq(spec)].iloc[0]
        a=infer_changes((q.post-q.pre).to_numpy(),q.oecd.to_numpy(),pre_months=int(m.pre_months),post_months=int(m.post_months),bootstrap=False)
        rows.append(dict(spec_id=spec,estimate=a['estimate'],CR1_SE=a['CR1_SE'],CR1_lower95=a['CR1_lower95'],CR1_upper95=a['CR1_upper95'],CR2_SE=a['CR2_SE'],CR2_lower95=a['CR2_lower95'],CR2_upper95=a['CR2_upper95'],CR2_p=a['CR2_p_two_sided']))
    actual=pd.DataFrame(rows);expected=pd.read_csv(ROOT/'expected/alternative_specification_results.csv')
    z=actual.merge(expected,on='spec_id',suffixes=('_actual','_expected'),validate='one_to_one');assert len(z)==2688
    errors={}
    for c in actual.columns[1:]:
        errors[c]=float(np.max(np.abs(z[c+'_actual']-z[c+'_expected'])))
        assert np.allclose(z[c+'_actual'],z[c+'_expected'],atol=1e-9,rtol=1e-7),c
    actual.to_csv(TABLES/'all_alternative_specification_results.csv',index=False)
    linked=actual.merge(definitions,on='spec_id');u=linked.drop_duplicates('effective_design_id')
    original=u[u.outcome.eq('FDAage2023median361')]
    summary=dict(labelled_specifications=len(actual),distinct_specifications=len(u),original_age_family_distinct=len(original),original_age_family_negative=int(original.estimate.lt(0).sum()),all_distinct_negative_fraction=float(u.estimate.lt(0).mean()))
    assert summary['distinct_specifications']==2100 and summary['original_age_family_distinct']==300 and summary['original_age_family_negative']==288
    print('Recomputed all 2,688 specifications (2,100 distinct).',flush=True)

    # Legacy model estimates and clustered covariance, from country-level
    # X'X and X'y totals. No medicine-level input is needed for this step.
    models=pd.read_csv(ROOT/'data/regression_model_definitions.csv')
    moments=pd.read_csv(ROOT/'data/regression_cluster_moments.csv.gz')
    result=[]
    for m in models.itertuples():
        q=moments[moments.model_id.eq(m.model_id)];k=int(m.parameters);AA=[];vv=[]
        for c,s in q.groupby('country',sort=False):
            A=np.zeros((k,k));v=np.zeros(k)
            for r in s.itertuples():
                i,j=int(r.term_i),int(r.term_j)
                if j==-1:v[i]=r.crossproduct
                else:A[i,j]=r.crossproduct;A[j,i]=r.crossproduct
            AA.append(A);vv.append(v)
        total=np.sum(AA,axis=0);rhs=np.sum(vv,axis=0);bread=np.linalg.pinv(total);b=bread@rhs
        score=[v-A@b for v,A in zip(vv,AA)];meat=sum(np.outer(v,v) for v in score)
        G=len(AA);cov=bread@meat@bread*G/(G-1)*(m.N-1)/(m.N-m.rank)
        j=int(m.target_column);se=float(np.sqrt(max(cov[j,j],0)));crit=t.ppf(.975,G-1)
        result.append(dict(model_id=m.model_id,family=m.family,label=m.label,estimate=b[j],cluster_SE=se,lower95=b[j]-crit*se,upper95=b[j]+crit*se))
    rr=pd.DataFrame(result);reference=pd.read_csv(ROOT/'expected/regression_results.csv');zz=rr.merge(reference,on='model_id',suffixes=('_actual','_expected'),validate='one_to_one')
    legacy_errors={}
    for c in ['estimate','cluster_SE','lower95','upper95']:
        legacy_errors[c]=float(np.max(np.abs(zz[c+'_actual']-zz[c+'_expected'])))
        assert np.allclose(zz[c+'_actual'],zz[c+'_expected'],atol=1e-7,rtol=1e-6),('regression moments',c,legacy_errors[c])
    rr.to_csv(TABLES/'main_and_historical_regressions.csv',index=False)
    primary=rr[rr.label.eq('Main six-country FDA-age balance')].iloc[0]
    assert abs(primary.estimate+.2610798097212942)<1e-9
    print('Recomputed',len(rr),'main and historical regressions from aggregate cross-products.',flush=True)

    cohort_input=pd.read_csv(ROOT/'data/Approval_cohort_country_outcomes.csv')
    cohort_reference=pd.read_csv(ROOT/'data/Approval_cohort_results.csv');cohort_rows=[]
    for cohort,q in cohort_input.groupby('cohort',sort=False):
        fit=infer_changes((q.post-q.pre).to_numpy(),q.country.isin(['Australia','Switzerland','United Kingdom']).to_numpy(),bootstrap=False)
        ref=cohort_reference[cohort_reference.cohort.eq(cohort)].iloc[0]
        for new,old in [('estimate','estimate'),('CR1_SE','cluster_SE'),('CR1_lower95','lower95'),('CR1_upper95','upper95')]:assert abs(fit[new]-ref[old])<1e-9,(cohort,new)
        cohort_rows.append(dict(cohort=cohort,estimate=fit['estimate'],cluster_SE=fit['CR1_SE'],lower95=fit['CR1_lower95'],upper95=fit['CR1_upper95']))
    pd.DataFrame(cohort_rows).to_csv(TABLES/'approval_cohort_results.csv',index=False)

    # Calendar checks use the published aggregate monthly figure table.
    monthly=pd.read_csv(ROOT/'data/calendar_country_month_panel.csv',parse_dates=['month']);trials=[]
    cc=['Australia','Switzerland','United Kingdom','Brazil','Hong Kong SAR','Singapore'];treat=np.array([1,1,1,0,0,0])
    windows=[('May2025_announcement',('2023-01-01','2025-04-01'),('2025-06-01','2026-05-01')),('May2024_placebo',('2023-01-01','2024-04-01'),('2024-06-01','2025-04-01')),('Jan2025_placebo',('2023-01-01','2024-12-01'),('2025-01-01','2025-04-01')),('Jan2026_vs_preannouncement',('2023-01-01','2025-04-01'),('2026-01-01','2026-05-01'))]
    for endpoint,q in monthly.groupby('endpoint'):
        for outcome in ['aggregate','exposure']:
            for label,(a,b),(c,d) in windows:
                pre=q[q.month.between(a,b)].groupby('country')[outcome].mean().reindex(cc)
                post=q[q.month.between(c,d)].groupby('country')[outcome].mean().reindex(cc)
                fit=infer_changes((post-pre).to_numpy(),treat,pre_months=len(pd.date_range(a,b,freq='MS')),post_months=len(pd.date_range(c,d,freq='MS')),bootstrap=True)
                trials.append(dict(test=label,endpoint=endpoint,outcome=outcome,**fit))
            slopes=[]
            for c in cc:
                z=q[q.country.eq(c)&q.month.between('2023-01-01','2025-04-01')].sort_values('month')
                slopes.append(np.polyfit(np.arange(len(z)),z[outcome].to_numpy(),1)[0])
            slopes=np.asarray(slopes);a=slopes[:3];b=slopes[3:];va=a.var(ddof=1)/3;vb=b.var(ddof=1)/3
            se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/2+vb**2/2);point=a.mean()-b.mean();crit=t.ppf(.975,df)
            trials.append(dict(test='Pre2025_linear_slope_difference',endpoint=endpoint,outcome=outcome,estimate=point,CR2_SE=se,CR2_lower95=point-crit*se,CR2_upper95=point+crit*se,WCR_webb_p_two_sided=np.nan))
    ct=pd.DataFrame(trials);ex=pd.read_csv(ROOT/'expected/calendar_tests.csv')
    z=ct.merge(ex,on=['test','endpoint','outcome'],suffixes=('_actual','_expected'),validate='one_to_one');assert len(z)==20
    for c in ['estimate','CR2_SE','CR2_lower95','CR2_upper95','WCR_webb_p_two_sided']:
        assert np.allclose(z[c+'_actual'],z[c+'_expected'],atol=1e-9,rtol=1e-7,equal_nan=True),('calendar',c)
    ct.to_csv(TABLES/'calendar_checks.csv',index=False)
    coverage=pd.read_csv(ROOT/'data/FDA_age_documented_funding_curves.csv')
    assert np.allclose(100*coverage.documented_by_horizon/coverage.denominator,coverage.documented_percent)
    coverage.to_csv(TABLES/'FDA_age_documented_funding_curves.csv',index=False)
    projection_checks=projections(ROOT,TABLES)
    proc=subprocess.run([sys.executable,str(ROOT/'code/make_figures.py')],cwd=ROOT,capture_output=True,text=True,encoding='utf8',errors='replace')
    (OUT/'figure_run.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf8');assert proc.returncode==0,proc.stderr
    report=dict(status='passed',paper='U.S. most-favored-nation pricing and foreign drug funding',SSRN='https://ssrn.com/abstract=7652061',doi='10.2139/ssrn.7652061',python=platform.python_version(),platform=platform.platform(),input_level='Final aggregate country-period, country-month, and regression cross-product tables; no source microdata',**summary,alternative_specification_maximum_errors=errors,main_and_historical_regressions=len(rr),approval_cohort_regressions=len(cohort_rows),regression_moment_maximum_errors=legacy_errors,calendar_comparisons=20,projection_tables=projection_checks,source_downloads_during_replication=False,source_archives_or_medicine_records_in_package=False)
    (OUT/'replication_report.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
