import xarray as xr
import numpy as np
import pandas as pd
###

import phippery
from phippery.normalize import svd_rank_reduction,enrichment, standardized_enrichment, differential_selection_wt_mut, rank_data, svd_aa_loc, size_factors, cpm
from phippery.collapse import collapse_sample_groups, pairwise_correlation_by_sample_group
from phippery.utils import *
#from phippery.tidy import tidy_ds
#from phippery.modeling import neg_binom_model


#import pickle
#from collections import defaultdict
#import copy
#from numpy.linalg import svd
#import scipy.stats as ss
#import matplotlib.pyplot as plt
#from matplotlib.colors import to_hex
#from matplotlib import cm
#import seaborn as sns
#import time
#import warnings
#from adjustText import adjust_text
import sys
#warnings.filterwarnings("ignore")

ds = phippery.load(sys.argv[1])

# a few things to fix up in the sample table
s1_s2 = id_coordinate_subset(ds, table="peptide_table", where="Protein", is_in=["S1","S2"])
ds = ds.loc[dict(sample_id=list(sams), peptide_id=s1_s2)]

pep_df = ds.peptide_table.to_pandas()
e = max(pep_df[pep_df["Protein"]=="S1"]["Loc"])
pep_df.loc[pep_df["Protein"]=="S2", "Loc"] += e
ds["peptide_table"] = xr.DataArray(pep_df, dims=ds.peptide_table.dims)

bat_ds = []
for batch, batch_ds in iter_sample_groups(ds, "library_batch"):
    
    batch_beads = id_coordinate_subset(batch_ds, where="control_status", is_equal_to="beads_only")
    batch_lib = id_coordinate_subset(batch_ds, where="control_status", is_equal_to="library")
    emp = id_coordinate_subset(batch_ds, where="control_status", is_equal_to="empirical")
    
    ret_ds = cpm(batch_ds, per_sample=True, data_table="counts", new_table_name="cpm", inplace=False)
    enrichment(ret_ds, list(batch_lib), "counts", new_table_name="counts_enrichment")
    standardized_enrichment(ret_ds, list(batch_lib), list(batch_beads), "counts", new_table_name="counts_std_enrichment")
    emp_beads = set.union(set(batch_beads), set(emp))

    emp_beads_ds = ret_ds.loc[dict(sample_id=list(emp_beads))]
    size_factors(emp_beads_ds)

    #beads_ds = emp_beads_ds.loc[dict(sample_id=batch_beads)]
    emp_ds = emp_beads_ds.loc[dict(sample_id=emp)]

    #neg_binom_model(
    #    emp_ds,
    #    beads_ds,
    #    nb_p=2,
    #    trim_percentile=99.9,
    #    data_table="size_factors",
    #    inplace=True,
    #    new_table_name="neg_binom_mlxp"
    #)

    differential_selection_wt_mut(
            emp_ds, 
            data_table="counts_enrichment", 
            relu_bias=None, 
            scaled_by_wt=True, 
            smoothing_flank_size=1,
            new_table_name="smooth_flank_1_enr_diff_sel"
    )

    bat_ds.append(emp_ds)
    
merged = bat_ds[0].merge(bat_ds[1])
phippery.dump(merged, f"{sys.argv[2]}")
