


import phippery
from phippery.utils import *
from phippery.tidy import tidy_ds

# load the dataset post stats computation
ds = phippery.load("../_ignore/curated-samples-06-29-21-bowtie2-layered.phip")
st = ds.sample_table.to_pandas()
st.loc[st[st["visit_number"]=="9"].index, "visit_number"] = "36 Days post-vaccination"
st.loc[st[st["visit_number"]=="12"].index, "visit_number"] = "119 Days post-vaccination"
ds["sample_table"] = xr.DataArray(st, dims=ds.sample_table.dims)

pt = ds.peptide_table.to_pandas()
site_wt = {row["Loc"]: row["aa_sub"] for idx, row in pt.iterrows() if row["is_wt"]==True}
#print(site_wt)
#import sys
#sys.exit()

condition_dfs = []
condition_groupby = ["participant_ID", "library_batch", "visit_number"]
#count = 0
for condition, condition_ds in phippery.iter_sample_groups(ds, condition_groupby):
    ds_slim = condition_ds.loc[dict(sample_metadata=["sample_group"])]
    slim_tall = tidy_ds(ds_slim)

#Index(['peptide_id', 'sample_id', 'counts_enrichment', 'counts_std_enrichment',
#           'size_factors', 'smooth_flank_1_enr_diff_sel', 'cpm', 'counts', 'Oligo',
#                  'Virus', 'Protein', 'Loc', 'aa_sub', 'is_wt', 'sample_group'],
#                        dtype='object')

    wt = [site_wt[site] for site in slim_tall["Loc"]]
    pro = ["A B C" for i in range(len(slim_tall))]
    cond = [", ".join(condition) for i in range(len(slim_tall))] 
    condition_dfs.append(
        pd.DataFrame({
            'site': slim_tall["Loc"],
            'protein_site': slim_tall["Loc"],
            'label_site': slim_tall["Loc"],
            'site_WT_enrichment': slim_tall["counts_enrichment"],
            'mut_scaled_diff_sel': slim_tall["smooth_flank_1_enr_diff_sel"],
            'mutation': slim_tall["aa_sub"],
            'wildtype': wt,
            'protein': pro,
            'condition': cond
        })
    )
    #print(sample_df)
    #count += 1
    #if count > 3: break

from functools import reduce
merged_samples_df = reduce(
    lambda l, r: l.append(r, ignore_index=True, verify_integrity=False),
    condition_dfs,
)
merged_samples_df.to_csv("dms-view-data-table.csv", index=False)
#print(merged_samples_df)









