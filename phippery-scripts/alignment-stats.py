
from plotnine import *
import matplotlib.pyplot as plt
from phippery.utils import *
import seaborn as sns
import xarray as xr
import pickle
import sys

# this is a messy way to access the data from my pickle dump binary
# phippery includes some niceer wrapper functions for this stuff
# but we'll forget about that for now

###
# load the pickled binary dataset
###
ds = pickle.load(open(sys.argv[1], "rb"))
not_beads = id_coordinate_subset(ds, where="control_status", is_not_equal_to="beads_only")
ds = ds.loc[dict(sample_id=not_beads)]

# The Dataset is fairly simple. a consize explanation ca be found
# in the README here: https://github.com/matsengrp/phippery
# for now, we can just look at the sample table
# Unfortunately, my xarray approach is slightly flawed in that it a data variable
# can only be of a single datatype, and this means all my metadata is stored as "Objects"
# should really find a better solution at some point. the "infer_objects" is my hack

###
# from the xarray Dataset Object, grab the sample table
###
sample_table = ds.sample_table.to_pandas().infer_objects()

###
# alignment stats to visualize
###

# I want to start by visualizing a few different alignment stats in the following columns
# these all get added during the pipeline run

# "raw_total_sequenced:"
# "reads mapped:"
# "error rate:"
# "average quality:"
sample_table["coverage:"] = sample_table["reads mapped:"] / len(ds.peptide_id)
align_stats = [c for c in sample_table.columns if c.endswith(":")]

# generally, we want to facet the subplots by alignment stat
# and then one of a few helpful ways to split sample groups
# 1. "control_status"
# 2. "library_batch"
# 3. "seq_dir"
# 4. "experiment"

# an example of control_status
facet = "control_status"
n_sample_types = len(set(sample_table[facet]))
n_alignment_stats = len(align_stats)
log_scale = False

# seaborn <- preferred
fig, ax = plt.subplots(
        n_sample_types, 
        n_alignment_stats,
        figsize=[8,8]
)

for subplot_row, (facet, facet_df) in enumerate(sample_table.groupby(facet)):
    for subplot_col, align_stat in enumerate(align_stats):
        ax_s = ax[subplot_row, subplot_col]
        sns.histplot(
                data=facet_df, 
                x=align_stat, 
                log_scale=log_scale,
                ax = ax_s
        )
        t = "\n".join(align_stat[:-1].split(" "))
        if subplot_col == 0: ax_s.set_ylabel(f"{facet}\nsample counts")
        if subplot_col != 0: ax_s.set_ylabel("")
        if subplot_row == 0: ax_s.set_title(t)
        if subplot_row != n_sample_types-1: ax_s.set_xlabel("")
        ax_s.set_xlabel("")
        ax_s.set_xticklabels(ax_s.get_xticks(), rotation = 45)

fig.subplots_adjust(hspace=0.5, wspace=0.7)
fig.savefig(sys.argv[2])
