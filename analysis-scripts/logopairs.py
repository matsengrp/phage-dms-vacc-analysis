
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rcParams["font.family"] = "Helvetica"
import seaborn as sns
from statannot import add_stat_annotation
import logomaker

from collections import defaultdict
import argparse
import os
import sys

import phippery
from phippery.tidy import tidy_ds
from phippery import *
from epitopes import EPITOPES
from epitopes import PEPTIDE_FLANK

parser = argparse.ArgumentParser(description='')
parser.add_argument('-subgroup', type=str)
parser.add_argument('-dataset', type=str)
parser.add_argument('-out', type=str)
parser.add_argument('-batch', type=str)
args = parser.parse_args()

batch = args.batch
flank = "flank_1"
group=args.subgroup

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
    sublabel = "B"

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
else:

    cmap=plt.get_cmap("Set3")
    infect_vacc_colors = [cmap.colors[8], cmap.colors[9]]
    figsize=[10, 10]
    t1 = "Pre-vaccine draw" 
    t2 = "Post-vaccine draw 1"
    l1 = "Pre-vaccine draw \n(post infection)" 
    l2 = "Post-vaccine draw"
    sublabel = "A"

    # Group 1
    hos = set(id_coordinate_subset(
        ds, where="sample_group", 
        is_equal_to="Hospitalized Serum"
    ))
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
    group3 = set.intersection(hos, pre)


    all_comparison_samples = set.union(group1, group2, group3)
    batch_samples = set(
            id_coordinate_subset(ds, where="library_batch", is_equal_to=f"{batch}")
    )
    batch_comparison_samples = set.intersection(batch_samples, all_comparison_samples)

ds = ds.loc[dict(sample_id=list(batch_comparison_samples))]
wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds_wt = ds.loc[dict(peptide_id=wt)]
enrichments = ds_wt["counts_enrichment"].to_pandas()
s_table = ds.sample_table.to_pandas()
p_table = ds_wt.peptide_table.to_pandas()
p_table_m = ds.peptide_table.to_pandas()

for epitope, metadata in EPITOPES.items():

    # CTD region too big to plot like this.
    if epitope=="CTD": continue
    limits = metadata["limits"]
    bt = metadata["binding threshold"]

    # compute the wildtype binding
    pep_id = p_table.loc[p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]
    s_table[f"{epitope}_raw_agg"] = epitope_enrichment.sum(axis=0).values

    # cut out samples below threshold and non-epitope specific peptides
    samples_above_thresh = s_table[s_table[f"{epitope}_raw_agg"]>=bt].index.values
    pep_id_mut = p_table_m.loc[p_table_m["Loc"].isin(range(limits[0], limits[1])),:].index
    ds_cur = ds.loc[dict(sample_id=samples_above_thresh, peptide_id=pep_id_mut)]

    # counts the number of samples in each category
    n_g1 = len(id_coordinate_subset(
        ds_cur, where="visit_number", is_equal_to=t1
    ))
    n_g2 = len(id_coordinate_subset(
        ds_cur, where="visit_number", is_equal_to=t2
        ))

    # counts the number of participants with samples in both categories
    paired_participants = []
    for p, p_df in iter_sample_groups(ds_cur, "participant_ID"):
        if len(p_df.sample_id) == 2: paired_participants.append(p)
    paired_participants = paired_participants[:10]
    num_pairs = len(paired_participants)
    if num_pairs==0:continue

    width = 10
    mosaic = np.zeros([num_pairs+1, width])
    for i in range(num_pairs):
        mosaic[i, 0:5] = int(i*2)
        mosaic[i, 5:11] = int(i*2+1)
    mosaic[num_pairs, :] = num_pairs*2

    chem_legend_ax = num_pairs*2

    mosaic = mosaic.astype(int)

    # set up the mosaic
    figsize = [16, max(num_pairs, 5)]
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    axd = fig.subplot_mosaic(
        mosaic,
        gridspec_kw={
            "wspace": 1.2,
            "hspace": 0.25,
        },
        empty_sentinel=-1,
    )

    
    kw = dict(ha="center", va="center", fontsize=19, color="black")
    sublabel = "B" if epitope == "CTD-N" else sublabel 
    sublabel = "A" if epitope == "NTD" else sublabel 
    sublabel = "A" if group != "moderna" else sublabel
    a = mosaic[0,0]
    axd[a].text(-0.16, 1.10, f"{sublabel})", transform=axd[a].transAxes, **kw)
    #axd[a].text(

    axd[mosaic[0,0]].set_title(l1, size=18)
    axd[mosaic[0,5]].set_title(l2, size=18)

    paired_samples = id_coordinate_subset(
            ds, 
            where="participant_ID", 
            is_in=paired_participants
    )
    ds_paired = ds_cur.loc[dict(sample_id=paired_samples)]
    for i, (pid, pid_ds) in enumerate(iter_sample_groups(ds_paired, "participant_ID")):
        #continue

        pid_tall = tidy_ds(pid_ds)
        p1 = pid_tall[pid_tall["visit_number"]==t1]
        p2 = pid_tall[pid_tall["visit_number"]==t2]
        
        p1_wide = p1.pivot("Loc", "aa_sub", f"smooth_{flank}_enr_diff_sel")
        p2_wide = p2.pivot("Loc", "aa_sub", f"smooth_{flank}_enr_diff_sel")

        crp_logo = logomaker.Logo(p1_wide,
                  font_name='Mono',
                  ax=axd[mosaic[i, 0]],
                  flip_below=False)


        crp_logo = logomaker.Logo(p2_wide,
                  font_name='Arial Rounded MT Bold',
                  ax=axd[mosaic[i, 5]],
                  flip_below=False)

        axd[mosaic[i][0]].set_frame_on(False)
        axd[mosaic[i][5]].set_frame_on(False)
        
        if i == num_pairs - 1:
            for s in [0, 5]:
                axd[mosaic[i][s]].set_xticks(range(limits[0], limits[1]))
                for tick in axd[mosaic[i][s]].get_xticklabels():
                    tick.set_rotation(90)
        else:
            for s in [0, 5]:
                axd[mosaic[i][s]].get_xaxis().set_visible(False)

        
        axd[mosaic[i][0]].set_ylabel(f"{pid}")

        ax = axd[mosaic[i][0]]

        if group != "moderna":
            dpi = set(p1["days_post_infection"])
            assert len(dpi) == 1
            dpi = list(dpi)[0]
            tt1 = f'{dpi} Days post-infection'
            ax.text(0.85, 0.35, f'{tt1}', horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)

        for item in ([ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)
        ax.yaxis.label.set_fontsize(18)


        ax = axd[mosaic[i][5]]

        if group != "moderna":
            dpi = set(p2["days_post_infection"])
            dpv = set(p2["days_post_vaccination"])

            assert len(dpi) == 1
            dpi = list(dpi)[0]
            tt1 = f'{dpi} Days post-infection'

            assert len(dpv) == 1
            dpv = list(dpv)[0]
            tt2 = f'{dpv} Days post-vaccination'

            ax.text(0.85, 0.35, f'{tt1}\n{tt2}', horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)

        for item in ([ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)

    chemistry = {
        "red" : "Acidic",
        "blue" : "Basic",
        "black" : "Hydrophobic",
        "purple" : "Neutral",
        "green" : "Polar"
    }
    legend_elements = [
        patches.Patch(
             facecolor=col,
             edgecolor="black",
             label= chemistry[col],
        )
        for col in chemistry
    ]
    if sublabel == "B":
        axd[chem_legend_ax].legend(
                handles=legend_elements,
                loc="center",
                frameon=False,
                ncol=5,
                prop={'size': 12},
                bbox_to_anchor=[0.5,0.15]
        )
    axd[chem_legend_ax].axis("off")

    group_title_map = {
        "moderna" : "Moderna Trial Cohort",
        "haarvi" : "HAARVI Cohort"
    }

    epitope_title_map = {
        "NTD" : "N-Terminal Domain (NTD)",
        "CTD-N" : "C-Terminal Domain (CTD-N)",
        "FP" : "Fusion Peptide (FP)",
        "stem helix-HR2" : "Stem Helix-HR2 (SH-H)"
    }
    figtitle = f"{group_title_map[group]}\n{epitope_title_map[epitope]}"
    fig.suptitle(figtitle, fontsize=20)

    plt.savefig(f"{epitope}-{args.out}", dpi=600)
