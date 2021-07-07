import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import to_hex
import matplotlib.patches as patches
from matplotlib import cm
import seaborn as sns
from statannot import add_stat_annotation

import sys
import copy
import warnings
warnings.filterwarnings("ignore")

import phippery
from phippery.utils import *



######
#KNOBS
######

batch="SPIKE2"

metric = "counts_enrichment"
highlight="sample_group"
np.random.seed(24)
mode_agg="wt_only"
meta=f"04-27-21/"
ag="sum"
colormap="Set3"
epitope_colormap="Set3"

ds = phippery.load(sys.argv[1])
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

epitope_order = [
    "NTD",
    "CTD",
    "FP", 
    "HR2",
]

epitope_limits = {
    "NTD" : [285, 305],
    "CTD" : [540, 690],
    "FP" : [800,833],
    "HR2" : [1130,1165],
}

cmap=plt.get_cmap(epitope_colormap)
epitope_colors = {
    "none":"black",
}
for i, (e, _) in enumerate(epitope_limits.items()):
    epitope_colors[e] = cmap.colors[i+8]

batch_beads = set(id_coordinate_subset(
        ds, 
        where="control_status", 
        is_equal_to="beads_only"
))
batch_emp = set(id_coordinate_subset(
        ds, 
        where="control_status", 
        is_equal_to="empirical"
))
batch_lib = set(id_coordinate_subset(
        ds, 
        where="control_status", 
        is_equal_to="library"
))

wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds = ds.loc[dict(peptide_id=wt)]

############################################
################ PREP DATA #################
############################################

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

for epitope, limits in epitope_limits.items():
    pep_id = bxr_p_table.loc[bxr_p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    bxr_p_table.loc[pep_id, "epitope"] = epitope
    epitope_enrichment = enrichments.loc[pep_id, :]
    if ag == "sum":
        bxr_s_table[f"{epitope}_raw_agg"] = epitope_enrichment.sum(axis=0).values
    else:
        bxr_s_table[f"{epitope}_raw_agg"] = epitope_enrichment.mean(axis=0).values        

# NIH - 3A
s_table_nih = bxr_s_table[bxr_s_table["source"]=="NIH"]

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

# we will make a figure for each epitope region

fig = plt.figure(constrained_layout=True, figsize=[7, 7])
axd = fig.subplot_mosaic(
    """
    AB
    CD
    EF
    GH
    """,
    gridspec_kw={
        "wspace": 0.0,
        "hspace": 0.15,
    },
)

subplot_keys = [
        ["A", "B"], 
        ["C", "D"], 
        ["E", "F"],
        ["G", "H"]
]

for subplot_ins, epitope in zip(subplot_keys, epitope_order):
    

    hue_color = [min(1, i*0.6) for i in epitope_colors[epitope]]

    ##############################
    #F#
    ##############################

    a = sns.boxplot(
            data=s_table_con,
            x = "dpi_range", 
            y = f"{epitope}_raw_agg", 
            ax = axd[subplot_ins[0]],
            #color=cmap.colors[2],
            #color=epitope_colors[epitope],
            #color=cmap.colors[4],
            color='white',
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_con,
            x = "dpi_range", 
            y = f"{epitope}_raw_agg", 
            ax = axd[subplot_ins[0]],
            color=cmap.colors[4],
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.legend([],[],frameon=False)
    a.set_xlabel("Time post-infection (days)")
    ylabel = f"{epitope} Epitope \nS{limits[0]} - S{limits[1]} \n summed enrichment"
    a.set_ylabel(ylabel, size=9, fontdict={'color':hue_color})

    add_stat_annotation(
            a, 
            data=s_table_con, 
            x="dpi_range", 
            y=f"{epitope}_raw_agg", 
            box_pairs=[
                    ("0-60", "60-180"), 
                    ("0-60", "180-370"), 
                    ("60-180", "180-370")
            ],
            test='Mann-Whitney', 
            text_format='star', 
            comparisons_correction=None,
            loc='inside', 
            verbose=0
    )

    ##############################
    #G#
    ##############################

    s_table_post_vacc.loc[s_table_post_vacc["Vaccine"]=="Pfizer", "Vaccine"] = "BNT162b2"
    s_table_post_vacc.loc[s_table_post_vacc["Vaccine"]=="Moderna", "Vaccine"] = "mRNA-1273"

    a = sns.boxplot(
            data=s_table_post_vacc,
            x = "Vaccine", 
            y = f"{epitope}_raw_agg", 
            ax = axd[subplot_ins[1]],
            color="white",
            linewidth = 2.5
    )
    a = sns.swarmplot(
            data=s_table_post_vacc,
            x = "Vaccine", 
            y = f"{epitope}_raw_agg", 
            ax = axd[subplot_ins[1]], 
            hue="sample_group",
            palette=[cmap.colors[1], cmap.colors[2]]
    )
    a.tick_params(axis='x', labelrotation = 45)
    a.legend([],[],frameon=False)
    a.set_xlabel("Vaccine")
    ylabel = ""
    a.set_ylabel(ylabel)

    add_stat_annotation(
            a, 
            data=s_table_post_vacc, 
            x="Vaccine", 
            y=f"{epitope}_raw_agg", 
            box_pairs=[
                    ("BNT162b2", "mRNA-1273"), 
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

legend_elements = [
    patches.Patch(
         facecolor=cmap.colors[1],
         edgecolor="black", 
         label="No prior\ninfection"
    ),
    patches.Patch(
         facecolor=cmap.colors[2],
         edgecolor="black", 
         label="Prior\ninfection"
    )
]

axd["B"].legend(
        handles=legend_elements, 
        loc="upper right",
        bbox_to_anchor=[1.50, 1],
        frameon=False
)

fontsize = 15
kw = dict(ha="center", va="center", fontsize=fontsize, color="black")
axd["B"].set_title("Post-vaccination only")
axd["A"].set_title("Infected sample only")
axd["A"].text(-0.1, 1.2, f"A)", transform=axd["A"].transAxes, **kw)
axd["B"].text(-0.1, 1.2, f"B)", transform=axd["B"].transAxes, **kw)
for ax in ["A", "B", "C", "D", "E", "F"]:
    axd[ax].get_xaxis().set_visible(False)

plt.tight_layout()
fig.savefig(sys.argv[2])




