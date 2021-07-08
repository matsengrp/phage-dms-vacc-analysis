import phippery
from phippery.tidy import tidy_ds
from phippery.normalize import enrichment
from phippery.utils import *

import xarray as xr
import numpy as np

from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

import os
import sys
import argparse

from epitopes import EPITOPES

parser = argparse.ArgumentParser(description='')
parser.add_argument('-subgroup', type=str)
parser.add_argument('-dataset', type=str)
parser.add_argument('-out', type=str)
args = parser.parse_args()

batch = "SPIKE2"
outdir = f"../sandbox/wt-binding-dist"
figure_size = [10, 5]
group=args.subgroup

#if not os.path.exists(outdir): os.mkdir(outdir)
ds = phippery.load(args.dataset)

# set up the sample groups and some plotting vars
if group == "moderna":

    cmap=plt.get_cmap("Set3")
    infect_vacc_colors = [cmap.colors[4], cmap.colors[2]]
    figsize=[10, 16]
    t1 = "36 days post-vaccination"
    t2 = "119 days post-vaccination"
    l1 = "36 days post-vaccination"
    l2 = "119 days post-vaccination"

    st = ds.sample_table.to_pandas()
    st.loc[st[st["visit_number"]=="9"].index, "visit_number"] = t1
    st.loc[st[st["visit_number"]=="12"].index, "visit_number"] = t2
    ds["sample_table"] = xr.DataArray(st, dims=ds.sample_table.dims)

    # Group 1
    d_36 = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to=t1
    ))
    group1 = set(d_36)

    d_119 = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to=t2
    ))
    group2 = set(d_119)

    all_comparison_samples = set.union(group1, group2)
    batch_samples = set(
            id_coordinate_subset(ds, where="library_batch", is_equal_to=f"{batch}")
    )
    batch_comparison_samples = set.intersection(batch_samples, all_comparison_samples)

#"haarvi"
elif group=="haarvi-non-hos":

    cmap=plt.get_cmap("Set3")
    infect_vacc_colors = [cmap.colors[8], cmap.colors[9]]
    figsize=[10, 8]
    t1 = "Pre-vaccine draw" 
    t2 = "Post-vaccine draw 1"
    l1 = "Pre-vaccine draw\n(Non-hospitalized prior infection)" 
    l2 = "Post-vaccine draw 1"

    # Group 1
    #hos = set(id_coordinate_subset(
    #    ds, where="sample_group", 
    #    is_equal_to="Hospitalized Serum"
    #))
    non_hos = set(id_coordinate_subset(
        ds, where="sample_group", 
        is_equal_to="Non-Hospitalized Serum"
    ))
    pre = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to="Pre-vaccine draw"
    ))
    group1 = set.intersection(non_hos, pre)


    vacc_prior = set(id_coordinate_subset(
        ds, where="sample_group",
        is_equal_to="Vaccinated - (Prior infection)"
    ))
    post = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to="Post-vaccine draw 1"
    ))
    group2 = set.intersection(vacc_prior, post)


    all_comparison_samples = set.union(group1, group2)
    batch_samples = set(
            id_coordinate_subset(ds, where="library_batch", is_equal_to=f"{batch}")
    )
    batch_comparison_samples = set.intersection(batch_samples, all_comparison_samples)

elif group=="haarvi-hos":

    cmap=plt.get_cmap("Set3")
    infect_vacc_colors = [cmap.colors[8], cmap.colors[9]]
    figsize=[10, 8]
    t1 = "Pre-vaccine draw" 
    t2 = "Post-vaccine draw 1"
    l1 = "Pre-vaccine draw\n(Hospitalized prior infection)" 
    l2 = "Post-vaccine draw 1"

    hos = set(id_coordinate_subset(
        ds, where="sample_group", 
        is_equal_to="Hospitalized Serum"
    ))
    pre = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to="Pre-vaccine draw"
    ))
    group1 = set.intersection(hos, pre)


    vacc_prior = set(id_coordinate_subset(
        ds, where="sample_group",
        is_equal_to="Vaccinated - (Prior infection)"
    ))
    post = set(id_coordinate_subset(
        ds, where="visit_number",
        is_equal_to="Post-vaccine draw 1"
    ))
    group2 = set.intersection(vacc_prior, post)


    group1_ds = ds.loc[dict(sample_id=list(group1))]
    group1_par = get_all_sample_metadata_factors(group1_ds, "participant_ID")
    group2_ds = ds.loc[dict(sample_id=list(group2))]
    group2_par = get_all_sample_metadata_factors(group2_ds, "participant_ID")
    all_parts = set.intersection(set(group1_par), set(group2_par))
    all_comparison_samples = id_coordinate_subset(
            ds, 
            where="participant_ID", 
            is_in=list(all_parts)
    )

    #all_comparison_samples = set.union(group1, group2)
    batch_samples = set(
            id_coordinate_subset(ds, where="library_batch", is_equal_to=f"{batch}")
    )
    batch_comparison_samples = set.intersection(batch_samples, all_comparison_samples)

ds = ds.loc[dict(sample_id=list(batch_comparison_samples))]
wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds_wt = ds.loc[dict(peptide_id=wt)]


#def identify_axes(ax_dict, fontsize=48):
#    """
#    Helper to identify the Axes in the examples below.
#
#    Draws the label in a large font in the center of the Axes.
#
#    Parameters
#    ----------
#    ax_dict : dict[str, Axes]
#        Mapping between the title / label and the Axes.
#    fontsize : int, optional
#        How big the label should be.
#    """
#    kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
#    for k, ax in ax_dict.items():
#        ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)

# set up ploting mosaic
num_ep = len(EPITOPES) -1
mosaic = np.zeros([num_ep, 2])
for i in range(num_ep):
    mosaic[i, 0] = i*2
    mosaic[i, 1] = i*2+1

fig = plt.figure(constrained_layout=True, figsize=[7, 7])
axd = fig.subplot_mosaic(
    mosaic,
    gridspec_kw={
        "wspace": 0.25,
        "hspace": 0.25,
    },
)

# grab the necessary data for plotting
enrichments = ds_wt["counts_enrichment"].to_pandas()
s_table = ds_wt.sample_table.to_pandas()
p_table = ds_wt.peptide_table.to_pandas()

# rearrange some things for plotting sake
s_table["visit_number"] = s_table["visit_number"].astype("category")
s_table["visit_number"] = s_table["visit_number"].cat.reorder_categories(
        [t1, t2]
)


# we're going to plot each epitope on a row
i = 0
for epitope, metadata in EPITOPES.items():

    # CTD region too big to plot like this.
    if epitope=="CTD": continue
    limits = metadata["limits"]
    bt = metadata["binding threshold"]

    # subset to out epitope region
    pep_id = p_table.loc[p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]

    # sum the enrichments for all samples in the epitope region
    s_table[f"{epitope}_raw_agg"] = epitope_enrichment.sum(axis=0).values
    
    # plot each paried set next to each other in columns
    for j, (visit, visit_st) in enumerate(s_table.groupby("visit_number")):

        # sort and re-order the values
        sg_ds_s = visit_st.sort_values(f"{epitope}_raw_agg", ascending=False).reset_index()
        order = sg_ds_s["participant_ID"].values
        sg_ds_s["participant_ID"] = sg_ds_s["participant_ID"].astype("category")
        sg_ds_s["participant_ID"] = sg_ds_s["participant_ID"].cat.reorder_categories(order)

        # plot each wiltype epitope binding
        sns.barplot(
                x="participant_ID", 
                y=f"{epitope}_raw_agg", 
                data=sg_ds_s,
                ax=axd[mosaic[i, j]],
                color="slategrey"
        )
        axd[mosaic[i, j]].axhline(bt)
        axd[mosaic[i, j]].get_xaxis().set_visible(False)
        ylabel = f"{epitope} Epitope \n{limits[0]}-{limits[1]}"
        if j==0: 
            axd[mosaic[i, j]].set_ylabel(ylabel)
        else: 
            axd[mosaic[i, j]].set_ylabel("")
    i+=1

axd[mosaic[0, 0]].set_title(l1)
axd[mosaic[0, 1]].set_title(l2)

for a in axd:
    ax = axd[a]
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(11)

kw = dict(ha="center", va="center", fontsize=11, color="black", rotation=90)
fig.text(0.02, 0.5, "Summed Enrichment", **kw)
plt.subplots_adjust(left=0.15)
#plt.tight_layout()
fig.savefig(args.out)
