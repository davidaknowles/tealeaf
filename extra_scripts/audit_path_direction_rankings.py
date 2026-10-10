#!/usr/bin/env python3
"""Cross existing EC local-path rankings with EC and junction directions (LR A100)."""
import json, numpy as np, pandas as pd
from pathlib import Path
from tealeaf.sc.replication_audit import ranked_direction_summary
D=Path('/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential')
t=pd.read_csv('analyses/path_direction_sources/variants/variant_directions.tsv.gz',sep='\t')
def shards(root):
    return pd.concat([pd.read_csv(p,sep='\t',usecols=['test_id','p_value','statistic']) for p in sorted(Path(root).glob('shard_*/paired_path.tsv'))])
ranks={'zeta64 calibrated':pd.read_csv(D/'local_path_full_testing_eb_pairwise/a64/paired_path.tsv',sep='\t',usecols=['test_id','p_value','raw_p_value','statistic']),
       'zeta64 raw':pd.read_csv(D/'local_path_full_testing_eb_pairwise/a64/paired_path.tsv',sep='\t',usecols=['test_id','raw_p_value','statistic']).rename(columns={'raw_p_value':'p_value'}),
       'zeta2 raw':pd.read_csv(D/'local_path_full_testing_eb_pairwise/a2/paired_path.tsv',sep='\t',usecols=['test_id','raw_p_value','statistic']).rename(columns={'raw_p_value':'p_value'}),
       'zeta1 fixed':shards(D/'local_path_pairwise_effect_usage_a1_20260903'),
       'zeta1 profiled':shards(D/'local_path_pairwise_usage_a1_profiled_20261010'),
       'zeta1 free':shards(D/'local_path_pairwise_usage_a1_free_20261010')}
dirs={'EC fixed':'ec_production_dot','EC free':'ec_free_dot','junction subject':'junction_subject_dot','junction pooled':'junction_pooled_dot'}
base=t.loc[t.lr_complete & t.lr_depth.ge(20)]
rows=[]
for universe,u in (('LR eligible, EC direction',base),('also junction identifiable',base.loc[np.isfinite(base.junction_pooled_dot)&base.junction_pooled_dot.ne(0)&np.isfinite(base.junction_subject_dot)&base.junction_subject_dot.ne(0)])):
    for rn,r in ranks.items():
        cols=['p_value']+(['raw_p_value'] if 'raw_p_value' in r else [])
        m=u.drop(columns=['p_value','raw_p_value','statistic']).merge(r,on='test_id')
        m=m.sort_values(cols+['statistic','test_id'],ascending=[True]*len(cols)+[False,True],kind='stable')
        for dn,dc in dirs.items():
            loc=m.loc[np.isfinite(m[dc])&m[dc].ne(0)].copy()
            if universe.startswith('LR') and dn.startswith('junction'): continue
            loc['rank']=np.arange(1,len(loc)+1); loc['method']=dn; loc['pooled_replicated']=loc[dc].gt(0)
            for s in ranked_direction_summary(loc,cutoffs=(100,)):
                rows.append({'universe':universe,'ranking':rn,'direction':dn,'n':s['n_available'],'A100':round(s['normalized_auc'],3),'c100':s['agreement']})
out=pd.DataFrame(rows); print(out.to_string(index=False)); out.to_csv('analyses/path_direction_sources/variants/rank_by_direction.tsv',sep='\t',index=False)
