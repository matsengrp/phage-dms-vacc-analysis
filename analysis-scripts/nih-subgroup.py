
import sys
import copy
import warnings
warnings.filterwarnings("ignore")
import argparse

import phippery
from phippery.utils import *
from epitopes import EPITOPES

import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation

######
#KNOBS
######

parser = argparse.ArgumentParser(description='')
parser.add_argument('-enrichment_metric', type=str, default="counts_enrichment")
parser.add_argument('-dataset', type=str)
parser.add_argument('-out', type=str)
parser.add_argument('-batch', type=str, default="SPIKE2")
args = parser.parse_args()

batch=args.batch
metric = args.enrichment_metric
show_all_p = False
highlight="sample_group"
ag="sum"
colormap="Set3"
epitope_colormap="Set3"

ds = phippery.load(args.dataset)
batch_samples = id_coordinate_subset(ds, where="library_batch", is_equal_to=batch)
ds = ds.loc[dict(sample_id=batch_samples)]

sample_group_order = [
        'Vaccinated - (No prior infection) Cohort 1', 
]

sample_group_labels = [
        'Vaccinated\n(No prior infection) - Cohort 1', 
]

cmap=plt.get_cmap(colormap)
sample_group_colors = dict(zip(sample_group_order, cmap.colors))

cmap=plt.get_cmap(epitope_colormap)
epitope_colors = {
    "none":"black",
}

i=0
for e, _ in EPITOPES.items():
    if e in ["CTD-N", "CTD-2", "CTD-3"]: continue
    epitope_colors[e] = cmap.colors[i+8]
    i += 1

wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds = ds.loc[dict(peptide_id=wt)]
enr_label = "Fold enrichment"

# Grab enrichments
enrichments = ds[f"{metric}"]

# Grab sample table
s_table = ds.sample_table.to_pandas()

# Grab peptide table
p_table = ds.peptide_table.to_pandas()
p_table["Loc"] = p_table["Loc"].astype(int) 
p_table["epitope"] = "none"

############################################
################ BOXPLOTS ##################
############################################

bxr_s_table = copy.deepcopy(s_table)
bxr_p_table = copy.deepcopy(p_table)

#for epitope, limits in epitope_limits.items():
for epitope, metadata in EPITOPES.items():
    if epitope in ["CTD-1", "CTD-2", "CTD-3"]: continue
    limits = metadata["limits"]
    pep_id = bxr_p_table.loc[bxr_p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    bxr_p_table.loc[pep_id, "epitope"] = epitope
    epitope_enrichment = enrichments.loc[pep_id, :]
    if ag == "sum":
        bxr_s_table[f"{epitope}_raw_agg"] = epitope_enrichment.sum(axis=0).values
    else:
        bxr_s_table[f"{epitope}_raw_agg"] = epitope_enrichment.mean(axis=0).values        


# NIH - 3A
s_table_nih = bxr_s_table[bxr_s_table["source"]=="NIH"]
s_table_nih.to_csv("nih_wt_enrichment_sum.csv", index=True)

s_table_nih.loc[s_table_nih["age_range"]=="[18-55) y.o.", "age_range"] = "18-55"
s_table_nih.loc[s_table_nih["age_range"]=="[55-70) y.o.", "age_range"] = "55-70"
s_table_nih.loc[s_table_nih["age_range"]=="[70-100) y.o.", "age_range"] = "70-100"
s_table_nih["age_range"] = s_table_nih["age_range"].astype("category")
s_table_nih["age_range"] = s_table_nih["age_range"].cat.reorder_categories(
        ["18-55", "55-70", "70-100"]
)
s_table_nih.loc[s_table_nih["dose"]=="100ug", "dose"] = "100"
s_table_nih.loc[s_table_nih["dose"]=="250ug", "dose"] = "250"
s_table_nih["dose"] = s_table_nih["dose"].astype("category")
s_table_nih["dose"] = s_table_nih["dose"].cat.reorder_categories(
        ["100", "250"]
)

s_table_nih["days_post_vaccination"] = s_table_nih["days_post_vaccination"].astype("category")
s_table_nih["days_post_vaccination"] = s_table_nih["days_post_vaccination"].cat.reorder_categories(
        ["36", "119"]
)

# Post vaccine draw
s_table_post_vacc = bxr_s_table[bxr_s_table["visit_number"]=="Post-vaccine draw 1"]

# Convalescent
s_table_con = bxr_s_table[
        (bxr_s_table["days_post_infection"]!="none") &
        (bxr_s_table["Vaccine"]=="none")
]
s_table_con["dpi"] = s_table_con["days_post_infection"].astype(int)
#print(len(s_table_con))
s_table_con["dpi_range"] = "0-60"
s_table_con.loc[s_table_con["dpi"]>=60, "dpi_range"] = "60-180"
s_table_con.loc[s_table_con["dpi"]>=180, "dpi_range"] = "180-370"
s_table_con["dpi_range"] = s_table_con["dpi_range"].astype("category")
s_table_con["dpi_range"] = s_table_con["dpi_range"].cat.reorder_categories(
        ["0-60", "60-180", "180-370"]
)


####################
# PLOTTING
####################

# Create the boxplots
fig = plt.figure(constrained_layout=True, figsize=[12, 8])
axd = fig.subplot_mosaic(
    """
    AAA.BBDD.CCEE
    KKK.LLNN.MMOO
    FFF.GGII.HHJJ
    PPP.QQSS.RRTT
    """,
    gridspec_kw={
        "wspace": 0.0,
        "hspace": 0.55,
    },
)

epitope_order = [
    "NTD",
    "CTD",
    "FP", 
    "stem helix-HR2",
]

subplot_keys = [
        ["A", "B", "C", "D", "E"],
        ["K", "L", "M", "N", "O"],
        ["F", "G", "H", "I", "J"],
        ["P", "Q", "R", "S", "T"]
]

for subplot_ins, epitope in zip(subplot_keys, epitope_order):
    

    ##############################
    #E#
    ##############################

    hue_color = [min(1, i*0.6) for i in epitope_colors[epitope]]
    a = sns.boxplot(
            data=s_table_nih, 
            x = "days_post_vaccination", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[0]], 
            palette = [epitope_colors[epitope], hue_color],
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_nih, 
            x = "days_post_vaccination", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[0]], 
            color=".25"
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.legend(bbox_to_anchor=(1.05, 1), loc=2, frameon=False)
    limits=EPITOPES[epitope]["limits"]
    a.set_title(f"")
    ylabel = f"{epitope} Epitope \nS{limits[0]} - S{limits[1]} \n"
    a.set_ylabel(ylabel)
    a.set_xlabel("")

    if show_all_p or epitope in ["CTD", "Linker-HR2"]:
        add_stat_annotation(
                a, 
                data=s_table_nih, 
                x="days_post_vaccination", 
                y=f"{epitope}_raw_agg", 
                box_pairs=[
                    ("36", "119"), 
                ],
                test='Wilcoxon', 
                text_format='full', 
                comparisons_correction='bonferroni',
                loc='inside', 
                verbose=2
        )

    ##############################
    #A#
    ##############################
    ylim = [0, max(s_table_nih[f"{epitope}_raw_agg"])+20]

    a = sns.boxplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"], 
            x = "dose", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[1]], 
            color = epitope_colors[epitope],
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"], 
            x = "dose", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[1]], 
            color = ".25"
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.set_title("")
    a.set_xlabel("")
    a.set_ylabel("")
    a.set_ylim(ylim)

    if show_all_p:
        add_stat_annotation(
                a, 
                data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"], 
                x="dose", 
                y=f"{epitope}_raw_agg", 
                box_pairs=[
                        ("100", "250"), 
                ],
                test='Mann-Whitney', 
                text_format='star', 
                comparisons_correction='bonferroni',
                #stats_params={'num_comparisons':36},
                loc='inside', 
                verbose=0
        )

    ##############################
    #B#
    ##############################

    a = sns.boxplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"], 
            x = "dose", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[3]], 
            #color = cmap.colors[1],
            color = hue_color,
            #color = epitope_colors[epitope],
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"], 
            x = "dose", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[3]], 
            color = ".25"
    )
    a.tick_params(axis='x', labelrotation = 45)
    #limits=epitope_limits[epitope]
    a.set_title(f"")
    a.set_xlabel("Dose ($\mu$g)")
    a.set_ylabel("")
    a.legend(bbox_to_anchor=(1.05, 1), loc=2, frameon=False)
    a.set_ylim(ylim)

    if show_all_p:
        add_stat_annotation(
                a, 
                data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"], 
                x="dose", 
                y=f"{epitope}_raw_agg", 
                box_pairs=[
                        ("100", "250"), 
                ],
                test='Mann-Whitney', 
                text_format='star', 
                comparisons_correction='bonferroni',
                #stats_params={'num_comparisons':36},
                loc='inside', 
                verbose=0
        )
        

    ##############################
    #C#
    ##############################

    hue_color = [min(1, i*0.6) for i in epitope_colors[epitope]]
    a = sns.boxplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"] ,
            x = "age_range", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[2]], 
            color = epitope_colors[epitope],
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"] ,
            x = "age_range", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[2]], 
            color=".25",
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.legend(bbox_to_anchor=(1.05, 1), loc=2, frameon=False)
    a.set_xlabel("")
    a.set_ylabel("")
    a.set_ylim(ylim)

    if show_all_p:
        add_stat_annotation(
                a, 
                data=s_table_nih[s_table_nih["days_post_vaccination"]=="36"],
                x="age_range", 
                y=f"{epitope}_raw_agg", 
                box_pairs=[
                        ("18-55", "55-70"), 
                        ("18-55", "70-100"), 
                        ("55-70", "70-100")
                               ],
                test='Mann-Whitney', 
                text_format='star', 
                comparisons_correction=None,
                #stats_params={'num_comparisons':36},
                loc='inside', 
                verbose=0
        )


    ##############################
    #D#
    ##############################
    
    a = sns.boxplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"],
            x = "age_range", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[4]], 
            color = hue_color,
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"],
            x = "age_range", 
            y = f"{epitope}_raw_agg", 
            ax=axd[subplot_ins[4]], 
            color=".25",
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.legend(bbox_to_anchor=(1.05, 1), loc=2, frameon=False)
    a.set_xlabel("Patient age (years)")
    ylim = a.get_ylim()

    if show_all_p:
        add_stat_annotation(
                a, 
                data=s_table_nih[s_table_nih["days_post_vaccination"]=="119"],
                x="age_range", 
                y=f"{epitope}_raw_agg", 
                box_pairs=[
                        ("18-55", "55-70"), 
                        ("18-55", "70-100"), 
                        ("55-70", "70-100")
                ],
                test='Mann-Whitney', 
                text_format='star', 
                comparisons_correction=None,
                loc='inside', 
                verbose=0
        )

##############################
#Fig Additions#
##############################

fontsize = 15
kw = dict(ha="center", va="center", fontsize=fontsize, color="black")
axd["A"].text(-0.3, 1.2, f"A)", transform=axd["A"].transAxes, **kw)
axd["B"].text(-0.3, 1.2, f"B)", transform=axd["B"].transAxes, **kw)
axd["C"].text(-0.3, 1.2, f"C)", transform=axd["C"].transAxes, **kw)
#for ax in ["A", "B", "C", "D", "E", "G", "H", "K", "L", "M", "N", "O", "", "R"]:
for ax in ["A", "B", "D", "C", "E", "K", "L", "N", "M", "O", "F", "G", "I", "H", "J"]:
    axd[ax].get_xaxis().set_visible(False)
for ax in ["D", "N", "I", "S", "E", "O", "J", "T"]:
    axd[ax].get_yaxis().set_visible(False)
axd["Q"].set_xlabel("36 DPV")
axd["S"].set_xlabel("119 DPV")
axd["R"].set_xlabel("36 DPV")
axd["T"].set_xlabel("119 DPV")
fontsize = 14
kw = dict(ha="center", va="center", fontsize=fontsize, color="black")

fig.text(0.48, 0.95, "Vaccine Dose ($\mu$g)", **kw)
#fig.text(0.48, 0.94, "Vaccine Dose ($\mu$g)", **kw)
fig.text(0.78, 0.95, "Participant Age (years)", **kw)
fig.text(0.21, 0.95, "Days post-vaccination (DPV)", **kw)

kw = dict(ha="center", va="center", fontsize=fontsize, color="black", rotation=90)
fig.text(0.02, 0.5, "Summed Enrichment", **kw)

plt.tight_layout()
#fig.savefig(f"../figures/heatmap-boxplot/{metric}-{batch}-36d.pdf")

#fig.suptitle(f"{epitope} Epitope Region S{limits[0]} - S{limits[1]}")
fig.savefig(args.out)

