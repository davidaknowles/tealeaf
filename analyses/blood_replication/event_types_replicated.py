import pandas as pd, numpy as np
F='analyses/coverage_null_filter_sh/pbmc'
g=pd.read_csv(f'{F}/split_filtered/split_gene_tests.tsv.gz',sep='\t'); g=g[g.universe.eq('matched')]
parts=[]
for ch in pd.read_csv(f'{F}/comparator_tests_filtered.tsv.gz',sep='\t',chunksize=2_000_000):
    parts.append(ch[ch.method.isin(['rMATS (paired)','LeafCutter','Tealeaf'])])
r=pd.concat(parts); r.to_parquet('.cache/adhoc/blood_rep/tests.parquet')
print(r.groupby(['comparison','method','fold']).size())
for comp,meth in [('rMATS (paired)','rMATS (paired)'),('LeafCutter','LeafCutter')]:
    d=g[g.comparison.eq(comp)]; c=d[~d.method.eq('Tealeaf')].set_index('gene_id'); t=d[d.method.eq('Tealeaf')].set_index('gene_id')
    rc=set(c.index[c.conjunction_q_value<=.05]); rt=set(t.index[t.conjunction_q_value<=.05])
    x=r[r.method.eq(meth)&r.comparison.eq(comp)&r.gene_id.isin(rc)&r.fold.astype(str).isin(['0','1'])]
    best=x.loc[x.groupby(['gene_id','fold']).p_value.idxmin()].copy()
    best['type']=best.feature_id.astype(str).str.split(':').str[1] if meth.startswith('rMATS') else 'cluster'
    mono=best.pair_id.str.contains('monocyte|DC|Macroph',case=False)
    print(comp, 'replicated', len(rc), 'Tealeaf-also', len(rc&rt)); print(best.groupby('fold').type.value_counts().unstack(0))
    print(' best pair involves monocytes/DC:', round(mono.mean(),2))
    if meth.startswith('rMATS'):
        allt=r[r.method.eq(meth)&r.fold.astype(str).eq('0')].feature_id.astype(str).str.split(':').str[1].value_counts(normalize=True).round(3); print(' all rMATS fold-0 tests by type', allt.to_dict())
