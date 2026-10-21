"""Regenerate the paper's five figures from the adjacent frozen CSVs.

Run this file from any working directory. Requires pandas and matplotlib.
Figure 1 uses phase means; Figure B1 alone uses observed monthly histories.
"""
from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,dates as mdates
HERE=Path(__file__).resolve().parents[1]/'data'
FIGOUT=Path(__file__).resolve().parents[1]/'build/figures'
FIGOUT.mkdir(parents=True,exist_ok=True)
for filename in ['times.ttf','timesbd.ttf','timesi.ttf','timesbi.ttf']:
    p=Path('C:/Windows/Fonts')/filename
    if p.exists():font_manager.fontManager.addfont(str(p))
family='Times New Roman' if Path('C:/Windows/Fonts/times.ttf').exists() else 'DejaVu Serif'
plt.rcParams.update({'font.family':family,'font.size':9,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'axes.titlesize':10,'figure.titlesize':11,'legend.fontsize':8.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'path'})
counts=pd.read_csv(HERE/'main_country_counts.csv')
rates=pd.read_csv(HERE/'age_balanced_country_period_outcomes.csv')
results=pd.read_csv(HERE/'all_six_country_results.csv')
oecd={'Australia','Switzerland','United Kingdom'}
counts['group']=counts.country.map(lambda c:'OECD' if c in oecd else 'Non-OECD')
rates['group']=rates.country.map(lambda c:'OECD' if c in oecd else 'Non-OECD')
raw=counts.groupby('group')[['pre_monthly','post_monthly']].mean()
age=rates.groupby('group')[['pre','post']].mean()
colors={'OECD':'#24617f','Non-OECD':'#b2602e'}
def save(fig,name):
    for ext in ['pdf','png','svg']:fig.savefig(FIGOUT/(name+'.'+ext),dpi=220,bbox_inches=None)
    plt.close(fig)

# Figure 1 combines the headline and the exposure split without repeating panels.
split=pd.read_csv(HERE/'Exposure_split_group_period_rates.csv')
fig,axes=plt.subplots(1,3,figsize=(7.2,3.35))
for group,color,style in [('OECD',colors['OECD'],'-'),('Non-OECD comparison',colors['Non-OECD'],'--')]:
    key='OECD' if group=='OECD' else 'Non-OECD'
    values=[100,100*raw.loc[key,'post_monthly']/raw.loc[key,'pre_monthly']]
    axes[0].plot([0,1],values,color=color,ls=style,marker='o',ms=4,lw=1.7,label=group)
    axes[0].annotate(f'{values[1]:.0f}',(1,values[1]),xytext=(5,0),textcoords='offset points',va='center',color=color,fontsize=8.5)
axes[0].set_title('A. All matched\nmedicines',loc='left',fontweight='bold',fontsize=10)
axes[0].set_ylabel('Entry rate (before = 100)');axes[0].set_ylim(0,135);axes[0].set_xlim(-.14,1.45)
axes[0].axhline(100,color='#aaa',ls=':',lw=.7)
for ax,spending,title in zip(axes[1:],['High','Low'],['B. Higher U.S.\nspending','C. Lower U.S.\nspending']):
    for group,color,style in [('OECD',colors['OECD'],'-'),('Non-OECD comparison',colors['Non-OECD'],'--')]:
        z=split[split.group.eq(group)&split.spending_group.eq(spending)].set_index('phase')
        y=[z.loc['Pre','rate'],z.loc['Post','rate']]
        ax.plot([0,1],y,color=color,ls=style,marker='o',ms=4,lw=1.7)
        for x,value in enumerate(y):ax.annotate(f'{value:.3f}',(x,value),xytext=(0,7 if group=='OECD' else -12),textcoords='offset points',ha='center',color=color,fontsize=8)
    ax.set_title(title,loc='left',fontweight='bold',fontsize=10)
    ax.set_ylim(0,.46);ax.set_xlim(-.22,1.22)
axes[1].set_ylabel('Entries per 100 medicines monthly')
for ax in axes:
    ax.set_xticks([0,1],['Before','After']);ax.grid(axis='y',alpha=.12)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.02),ncol=2,frameon=False)
fig.suptitle('Funding entry before and after MFN',fontweight='bold',y=.98)
fig.subplots_adjust(left=.075,right=.98,bottom=.22,top=.75,wspace=.55)
save(fig,'Figure_1_Combined')

fig,ax=plt.subplots(figsize=(7.2,2.8))
names=['Main six-country FDA-age balance','Main six-country Uniform matched medicines','Main six-country Entropy-balanced drug characteristics']
labels=['Approval-age balance (main)','Equal medicine weights','Four-characteristic balance']
for i,name in enumerate(names):
    r=results[results.test.eq(name)].iloc[0]
    ax.plot([r.lower95,r.upper95],[i,i],color=colors['OECD'],lw=1.6)
    ax.scatter(r.estimate,i,color=colors['OECD'],s=21,zorder=3)
    ax.annotate(f'{r.estimate:+.3f}',(r.estimate,i),xytext=(0,7),textcoords='offset points',ha='center',fontsize=8.5)
ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set_ylim(2.45,-.55)
ax.axvline(0,color='#777',ls=':',lw=.8);ax.grid(axis='x',alpha=.13)
ax.set_xlabel('Relative OECD change in higher-minus-lower spending entries\n(per 100 medicines monthly)')
ax.set_title('The result under three weighting methods',fontweight='bold',pad=12,fontsize=11)
fig.subplots_adjust(left=.29,right=.98,bottom=.27,top=.80)
save(fig,'Figure_Alternative_Weights')

# Figure 2 answers the exploratory, predefined approval-cohort question.
cohorts=pd.read_csv(HERE/'Approval_cohort_results.csv')
fig,ax=plt.subplots(figsize=(7.2,3.0))
for i,r in enumerate(cohorts.itertuples()):
    ax.plot([r.lower95,r.upper95],[i,i],color=colors['OECD'],lw=1.7)
    ax.scatter(r.estimate,i,color=colors['OECD'],s=25,zorder=3)
    ax.annotate(f'{r.estimate:+.3f}',(r.estimate,i),xytext=(0,8),textcoords='offset points',ha='center',fontsize=9)
labels=[f'{r.cohort} ({r.High+r.Low} medicines)' for r in cohorts.itertuples()]
ax.set_yticks(range(len(cohorts)),labels);ax.invert_yaxis();ax.set_ylim(2.5,-.6)
ax.axvline(0,color='#777',ls=':',lw=.8);ax.grid(axis='x',alpha=.13)
ax.set_xlabel('Relative OECD change in higher-minus-lower spending entries\n(per 100 medicines within the FDA approval cohort monthly)')
ax.set_title('Relative funding changes by FDA approval cohort',fontweight='bold',fontsize=11,pad=14)
fig.subplots_adjust(left=.29,right=.98,bottom=.27,top=.80)
save(fig,'Figure_2_Approval_Cohorts')

# Accounting decomposition relative to the SAME comparison mean in every OECD market.
is_oecd=rates.country.isin(['Australia','Switzerland','United Kingdom'])
comparison_mean=rates.loc[~is_oecd,'change'].mean()
contrib=rates.loc[is_oecd,['country','change']].copy()
contrib['comparison_mean_change']=comparison_mean
contrib['relative_change']=contrib.change-comparison_mean
contrib['contribution_to_pooled_estimate']=contrib.relative_change/3
contrib=contrib.sort_values('contribution_to_pooled_estimate')
main=results[results.test.eq('Main six-country FDA-age balance')].iloc[0]
assert abs(contrib.contribution_to_pooled_estimate.sum()-main.estimate)<1e-12
contrib.to_csv(FIGOUT/'Country_contributions_main_estimate.csv',index=False)
donor=rates.loc[~is_oecd,['country','change']].copy().sort_values('country')
donor.to_csv(FIGOUT/'Comparison_country_changes.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(7.2,3.35))
for i,r in enumerate(contrib.itertuples()):
    value=r.contribution_to_pooled_estimate
    axes[0].barh(i,value,color=colors['OECD'],height=.5)
    axes[0].text(value-.004,i,f'{value:+.3f}',ha='right',va='center',fontsize=9)
axes[0].set_yticks(range(len(contrib)),contrib.country.replace({'United Kingdom':'Britain'}))
axes[0].invert_yaxis();axes[0].set_xlim(contrib.contribution_to_pooled_estimate.min()*1.34,.008)
axes[0].set_title('A. OECD contributions',loc='left',fontweight='bold')
axes[0].set_xlabel('Contribution to pooled estimate\n(entries per 100 medicines monthly)')
for i,r in enumerate(donor.itertuples()):
    axes[1].barh(i,r.change,color=colors['Non-OECD'],height=.5)
    axes[1].text(r.change+(.007 if r.change>=0 else -.007),i,f'{r.change:+.3f}',ha='left' if r.change>=0 else 'right',va='center',fontsize=9)
axes[1].set_yticks(range(len(donor)),donor.country.replace({'Hong Kong SAR':'Hong Kong'}))
axes[1].invert_yaxis();axes[1].set_xlim(-.09,.30)
axes[1].axvline(comparison_mean,color='#555',ls='--',lw=1,label=f'Comparison mean: {comparison_mean:+.3f}')
axes[1].set_title('B. Comparison-market changes',loc='left',fontweight='bold')
axes[1].set_xlabel('Change after announcement\n(entries per 100 medicines monthly)')
axes[1].legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(0,-.32))
for ax in axes:
    ax.axvline(0,color='#888',lw=.8);ax.grid(axis='x',alpha=.12)
fig.suptitle('Country contributions to the main estimate',fontweight='bold',y=.98)
fig.subplots_adjust(left=.13,right=.97,bottom=.27,top=.80,wspace=.64)
save(fig,'Figure_2_Country_Contributions')

history=pd.read_csv(HERE/'Synthetic_control_monthly_series.csv',parse_dates=['date']).set_index('date')
fig,ax=plt.subplots(figsize=(7.2,3.35))
smooth_records=[]
for column,label,color,style in [('OECD_age_balanced','OECD',colors['OECD'],'-'),('synthetic_age_balanced','Non-OECD comparison',colors['Non-OECD'],'--')]:
    # Full six-month windows within each regime. No partial windows, invented
    # months, cross-policy averaging or smoothing-parameter selection by result.
    for j,part in enumerate([history.loc[:'2025-04-01'],history.loc['2025-06-01':]]):
        smooth=part[column].rolling(6,min_periods=6).mean()
        ax.plot(smooth.index,smooth,color=color,ls=style,lw=2.1,label=label if j==0 else None)
        for date,value in smooth.dropna().items():smooth_records.append({'date':date,'series':column,'six_month_average':value,'phase':'Pre' if j==0 else 'Post'})
pd.DataFrame(smooth_records).to_csv(FIGOUT/'Synthetic_control_six_month_averages.csv',index=False)
ax.axvline(pd.Timestamp('2025-05-12'),color='#777',ls=':',lw=.9)
ax.annotate('May 2025\nMFN announcement',xy=(mdates.date2num(pd.Timestamp('2025-05-12')),.97),xycoords=ax.get_xaxis_transform(),xytext=(5,0),textcoords='offset points',ha='left',va='top',fontsize=9,color='#555')
ax.text(.015,.98,'Six-month moving averages',transform=ax.transAxes,va='top',fontsize=9,color='#555')
ax.axhline(0,color='#aaa',ls=':',lw=.7)
ax.set_xlim(pd.Timestamp('2023-01-01'),pd.Timestamp('2026-05-31'))
ax.xaxis.set_major_locator(mdates.YearLocator());ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.grid(axis='y',alpha=.13)
ax.set_ylabel('Higher-minus-lower U.S.-spending entries\n(per 100 medicines monthly)')
ax.set_title('Funding trends before and after MFN',fontsize=11,fontweight='bold',pad=14)
ax.legend(frameon=False,loc='upper left',ncol=1,bbox_to_anchor=(0,-.16))
fig.subplots_adjust(left=.11,right=.98,bottom=.27,top=.84)
save(fig,'Figure_Synthetic_OECD_Comparison')
(FIGOUT/'figure_font_settings.json').write_text(json.dumps({'font_family':family,'base_font_points':9,'canvas_width_inches':7.2,'figures':['Figure_1_Combined','Figure_2_Approval_Cohorts','Figure_2_Country_Contributions','Figure_Alternative_Weights','Figure_Synthetic_OECD_Comparison'],'Figure_1_footer_lines_removed':True,'monthly_synthetic_smoothing':'Fixed full six-month trailing average, separately pre and post; no May bridge'},indent=2),encoding='utf8')
