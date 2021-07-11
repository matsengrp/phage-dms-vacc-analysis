import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches
import seaborn as sns
from statannot import add_stat_annotation
#####

import phippery
from phippery.utils import *
from epitopes import EPITOPES

import sys
import warnings
warnings.filterwarnings("ignore")
import argparse

print(EPITOPES)

######
#KNOBS
######

parser = argparse.ArgumentParser(description='')
parser.add_argument('-enrichment_metric', type=str, default="counts_enrichment")
parser.add_argument('-dataset', type=str)
parser.add_argument('-out', type=str)
parser.add_argument('-batch', type=str, default="SPIKE2")
args = parser.parse_args()

metric = ars.enrichment_metric
highlight="sample_group"
batch = args.batch
np.random.seed(24)
mode_agg="wt_only"
ag="sum"
max_enrichment_value=150
sample_colormap="Set3"
epitope_colormap="Set1"
heatcolormap="RdPu"
window_size=10
pca_axes_p_table = "pca_p_table.csv"

ds = phippery.load(args.dataset)
batch_samples = id_coordinate_subset(ds, where="library_batch", is_equal_to=batch)
ds = ds.loc[dict(sample_id=batch_samples)]

hcmap = getattr(cm, heatcolormap)

sample_group_order = [
        'Vaccinated - (No prior infection) Cohort 1', 
        'Vaccinated - (No prior infection) Cohort 2', 
        'Vaccinated - (Prior infection) Cohort 2', 
        'Hospitalized Serum Cohort 2', 
        'Non-Hospitalized Serum Cohort 2', 
        'healthy control Cohort 2'
]

sample_group_labels = [
        'Vaccinated\n(No prior infection) - Moderna Trial', 
        'Vaccinated\n(No prior infection) - HAARVI Trial', 
        'Vaccinated\n(Prior infection)', 
        'Hospitalized Serum', 
        'Non-Hospitalized Serum', 
        'Naive Serum'
]

cmap=plt.get_cmap(sample_colormap)
sample_group_colors = dict(zip(sample_group_order, cmap.colors))

wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds = ds.loc[dict(peptide_id=wt)]
enr_label = "Enrichment"

############################################
################ CURATE SAMPLES ############
############################################

# remove second time point from nih moderna samples
time_point = set(id_coordinate_subset(
    ds, 
    where="days_post_vaccination", 
    is_equal_to="119"
))

# cohort
cohort = set(id_coordinate_subset(
    ds, 
    where="manuscript_cohort", 
    is_equal_to="Cohort 1"
))

cohort_119 = set.intersection(
        time_point,
        cohort
    )

sams = set(ds.sample_id.values)
sams = sams - cohort_119
ds = ds.loc[dict(sample_id=list(sams))]

############################################
################ PREP DATA #################
############################################

# Grab enrichments
enrichments = ds[f"{metric}"]

# Grab sample table
s_table = ds.sample_table.to_pandas()
s_table["sample_group"] = s_table["sample_group"] + " " + s_table["manuscript_cohort"]
s_table["sample_group"] = s_table["sample_group"].astype("category")
s_table["sample_group"] = s_table["sample_group"].cat.reorder_categories(sample_group_order)
s_table["vaccination_status"] = "vaccinated"
s_table.loc[s_table["Vaccine"]=="none","vaccination_status"] = "non-vaccinated"
s_table["vaccination_status"] = s_table["vaccination_status"].astype("category")
s_table["vaccination_status"] = s_table["vaccination_status"].cat.reorder_categories(
        ["vaccinated", "non-vaccinated"]
        )


# Grab peptide table
p_table = ds.peptide_table.to_pandas()
p_table["Loc"] = p_table["Loc"].astype(int) 
p_table["epitope"] = "none"

############################################
################ GRID LAYOUT ###############
############################################

fig = plt.figure(constrained_layout=True, figsize=[14, 10])
axd = fig.subplot_mosaic(
    """
    .....
    AAA..
    AAAHG
    AAABC
    AAABC
    AAABC
    FFFDE
    IIIDE
    JJJDE
    """,
    gridspec_kw={
        "wspace": 0.55,
        "hspace": 0.25,
    },
)

kw = dict(ha="center", va="center", fontsize=16, color="black")
for ax, sublabel in zip(["A", "F", "B"], ["A", "B", "C"]):
    axd[ax].text(-0.10, 1.10, f"{sublabel})", transform=axd[ax].transAxes, **kw)

############################################
################ BOXPLOTS ##################
############################################

for epitope, metadata in EPITOPES.items():
    if epitope in ["CTD-1", "CTD-2", "CTD-3"]: continue
    limits = metadata["limits"]
    pep_id = p_table.loc[p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    p_table.loc[pep_id, "epitope"] = epitope
    epitope_enrichment = enrichments.loc[pep_id, :]
    if ag == "sum":
        s_table[f"{epitope}_raw_agg"] = epitope_enrichment.sum(axis=0).values
    else:
        s_table[f"{epitope}_raw_agg"] = epitope_enrichment.mean(axis=0).values        

boxpairs = [
        ("Non-Hospitalized Serum Cohort 2",g) 
        for g in set(sample_group_order)-set(["Non-Hospitalized Serum Cohort 2"])
]

epitope_order = [
    "NTD",
    "CTD",
    "FP", 
    "HR2",
]

first=True
for ax, epitope in zip(["B", "C", "D", "E"], epitope_order):
    sns.boxplot(
            data=s_table, 
            x="sample_group", 
            y=f"{epitope}_raw_agg", 
            ax=axd[ax], 
            palette=sample_group_colors,
    )

    add_stat_annotation(
            axd[ax], 
            data=s_table, 
            x="sample_group", 
            y=f"{epitope}_raw_agg", 
            order=sample_group_order,
            box_pairs = boxpairs,
            test='Mann-Whitney', 
            comparisons_correction=None,
            text_format='star', 
            loc='inside', 
            verbose=0
    )
    limits = EPITOPES[epitope]["limits"]
    ylabel = "Summed enrichment"
    title = f"{epitope} Epitope \n{limits[0]} - {limits[1]}"
    axd[ax].set_ylabel(ylabel)
    axd[ax].set_xlabel("")
    axd[ax].set_title(title, fontdict={'color':'slategray'})
    axd[ax].get_xaxis().set_visible(False)
    first = False


############################################
################ HEATMAP ###################
############################################

windows = list(range(0, max(p_table["Loc"]), window_size))
enrichment_columns = []
for l in range(len(windows)-1):
    start = windows[l]
    end = windows[l+1]
    pep_id = p_table.loc[p_table["Loc"].isin(range(start, end)),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]
    enrichment_columns.append(f"[{start}, {end})")
    if ag == "sum":
        agg = epitope_enrichment.sum(axis=0).values
    else:
        agg = epitope_enrichment.mean(axis=0).values
    s_table[f"[{start}, {end})"] = [min(max_enrichment_value, v) for v in agg]


sample_order = []
for g, g_df in s_table.groupby(
        ["vaccination_status", "manuscript_cohort", "sample_group"]
        ):
    sample_order.extend(g_df.index)
    
s_table = s_table.reindex(sample_order)

s_e = s_table[enrichment_columns]

cbar_kws = dict(use_gridspec=False, location="right", label=enr_label)
sns.heatmap(s_e, ax = axd["A"], cmap=hcmap, cbar_kws=cbar_kws)
axd["A"].get_yaxis().set_visible(False)
axd["A"].set_ylabel("Amino acid position windows")


base=0
for g, g_df in s_table.groupby(["manuscript_cohort", "sample_group"]):
    
    height = len(g_df.index)
    label = None if g[0] == "Cohort 2" and g[1] == "Vaccinated - (No prior infection)" else g[1]


    rect_v = patches.Rectangle(
            (-3, base),
            width=3,
            height=height,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            label=label,
            facecolor=sample_group_colors[g[1]]
    )
    axd["A"].axhline(base + height, color="black")
    axd["A"].add_patch(rect_v)
    base = base + height


base=0
for g, g_df in s_table.groupby(["manuscript_cohort"]):
    
    height = len(g_df.index)

    if g == "Cohort 1":
        label = "Moderna Trial"
    else:
        label = "HAARVI Trial"

    rect_v = patches.Rectangle(
            (-10, base),
            width=7,
            height=height,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor="None"
    )
    axd["A"].axhline(base + height, lw=2, color="black")
    axd["A"].text(-6.0, (base+(base+height))/2, label, rotation=90, va="center", ha="center", size=14)
    axd["A"].add_patch(rect_v)
    #axd["A"].set_xlabel("Amino acid position")
    base = base + height

base=0
for g, g_df in s_table.groupby(["vaccination_status"]):
    
    height = len(g_df.index)

    rect_v = patches.Rectangle(
            (-17, base),
            width=7,
            height=height,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor="None"
    )
    if base == 0:
        axd["A"].axhline(base + height, lw=3, color="black", linestyle="--")
    axd["A"].text(-13.0, (base+(base+height))/2, g, rotation=90, va="center", ha="center", size=14)
    axd["A"].add_patch(rect_v)
    base = base + height

axd["A"].set_xlabel("Amino acid position")
protein_regions = {
    "NTD" : [29, 294],
    "RBD" : [332, 523],
    "CTD" : [530, 681],
    "FP" : [816, 833],
    "HR1" : [910, 987],
    "HR2" : [1162, 1204]
}


for pr, limits in protein_regions.items():

    start = limits[0]//10
    width = ((limits[1] - limits[0]) // 10) + 1
    middle = start + (width//2)
    anno = f"{pr}\n{limits[0]}- \n {limits[1]}"

    rect_v = patches.Rectangle(
            (start, -7),
            width=width,
            height=7,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor="white"
    )

    axd["A"].text(middle, -16, anno, va="center", ha="center", size=12)
    axd["A"].add_patch(rect_v)

for ep, metadata in EPITOPES.items():
    if epitope in ["CTD-1", "CTD-2", "CTD-3"]: continue

    limits = metadata["limits"]
    start = limits[0]//10
    width = ((limits[1] - limits[0]) // 10) + 1
    middle = start + (width//2)

    rect_v = patches.Rectangle(
            (start, -4),
            width=width,
            height=4,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor="slategrey"
    )

    axd["A"].add_patch(rect_v)

axd["A"].get_xaxis().set_visible(False)

legend_elements = [
    patches.Patch(
         facecolor=sample_group_colors[sg], 
         edgecolor="black", 
         label=sample_group_labels[i]
    )
    for i, sg in enumerate(sample_group_order)
]

axd["H"].axis('off')
axd["H"].legend(
        handles=legend_elements, 
        loc="upper right",
        frameon=False
)
axd["G"].axis('off')
axd["G"].text(0.95, 0.95,
        """
        ns: 5.00e-02 < p <= 1.00e+00
        *: 1.00e-02 < p <= 5.00e-02
        **: 1.00e-03 < p <= 1.00e-02
        ***: 1.00e-04 < p <= 1.00e-03
        ****: p <= 1.00e-04
        """, 
        size=10, rotation=0,
         ha="right", va="top",
         )

pca_p_table = pd.read_csv(pca_axes_p_table)
for i, lo in enumerate(["F", "I", "J"]):
    for e, metadata in EPITOPES.items():
        e_l = metadata["limits"]
        axd[lo].plot(
            pca_p_table["Loc"].astype(int).values[e_l[0]:e_l[1]],
            pca_p_table[f"loading-{i+1}"].values[e_l[0]:e_l[1]],
            c = "slategrey",
            lw = 4,
            label=f"{e} - [{e_l[0]}, {e_l[1]})"
        )
    axd[lo].plot(
        pca_p_table["Loc"].astype(int),
        pca_p_table[f"loading-{i+1}"],
        c = "lightgray",
        lw = 1.2
    )
    axd[lo].margins(x=0)

axd["F"].get_xaxis().set_visible(False)
axd["F"].yaxis.set_label_position("right")
axd["F"].set_ylabel("Comp 1")

axd["I"].get_xaxis().set_visible(False)
axd["I"].yaxis.set_label_position("right")
axd["I"].set_ylabel("Comp 2\nPrincipal axes/directions\nin feature space")

axd["J"].set_xlabel("Amino Acid Position")
axd["J"].yaxis.set_label_position("right")
axd["J"].set_ylabel("Comp 3")
axd["J"].yaxis.label.set_fontsize(10)

axd["C"].set_ylabel("")
axd["E"].set_ylabel("")

for a in axd:
    ax = axd[a]
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(14)

plt.tight_layout()
fig.savefig(args.out)
