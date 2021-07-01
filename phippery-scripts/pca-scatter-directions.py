
import sys
import warnings
warnings.filterwarnings("ignore")

import phippery
from phippery.utils import *
from epitopes import EPITOPES

import xarray as xr
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
from matplotlib import cm
import seaborn as sns


######
#KNOBS
######

metric = "counts_enrichment"
batch = "SPIKE2"
np.random.seed(24)
epitope_colormap="Set3"
sample_colormap="Set3"
pcs=3

ds = phippery.load(sys.argv[1])
batch_samples = id_coordinate_subset(ds, where="library_batch", is_equal_to=batch)
ds = ds.loc[dict(sample_id=batch_samples)]

sample_group_order = [
        'Vaccinated - (No prior infection) Cohort 1', 
        'Vaccinated - (No prior infection) Cohort 2', 
        'Vaccinated - (Prior infection) Cohort 2', 
        'Hospitalized Serum Cohort 2', 
        'Non-Hospitalized Serum Cohort 2', 
        #'healthy control Cohort 2'
]

sample_group_labels = [
        'Vaccinated\n(No prior infection) - Moderna Trial', 
        'Vaccinated\n(No prior infection) - HAARVI Trial', 
        'Vaccinated\n(Prior infection)', 
        'Hospitalized Serum', 
        'Non-Hospitalized Serum', 
]

cmap=plt.get_cmap(sample_colormap)
sample_group_colors = dict(zip(sample_group_order, cmap.colors))

cmap=plt.get_cmap(epitope_colormap)
epitope_colors = {
    "none":"black",
}

i=0
for e, _ in EPITOPES.items():
    if e in ["CTD-1", "CTD-2", "CTD-3"]: continue
    epitope_colors[e] = cmap.colors[i+8]
    i += 1

wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds = ds.loc[dict(peptide_id=wt)]

############################################
################ CURATE SAMPLES ############
############################################

# remove 250 dose samples
mg_250 = set(id_coordinate_subset(
    ds, 
    where="dose", 
    is_equal_to="250ug"
))

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
    #is_equal_to="Moderna Trial"
    is_equal_to="Cohort 1"
))

cohort_119 = set.intersection(
        mg_250,
        cohort
    )
cohort_250 = set.intersection(
        time_point,
        cohort
    )

sams = set(ds.sample_id.values)
sams = sams - set.union(cohort_119, cohort_250)
sams = sams - set(id_coordinate_subset(
    ds, 
    where="sample_group", 
    is_equal_to="healthy control"
))
sams = sams - set([1333, 244]) # outliers
outliers = set(id_coordinate_subset(ds, where="participant_ID", is_equal_to="163C"))
sams = sams - outliers
ds = ds.loc[dict(sample_id=list(sams))]

############################################
################ PREP DATA #################
############################################


# Grab enrichments
enrichments_pca = ds[f"{metric}"]
enrichments = ds[f"{metric}"]

# Grab sample table
s_table = ds.sample_table.to_pandas()
s_table["sample_group"] = s_table["sample_group"] + " " + s_table["manuscript_cohort"]
s_table["sample_group"] = s_table["sample_group"].astype("category")
s_table["sample_group"] = s_table["sample_group"].cat.reorder_categories(sample_group_order)

# Grab peptide table
p_table = ds.peptide_table.to_pandas()
p_table["Loc"] = p_table["Loc"].astype(int) 
p_table["epitope"] = "none"

# Compute Principal Componany Ananlysis via mean-centered SVD
############################################
################ PCA #######################
############################################

print(f"enrichments shape: {enrichments_pca}")
pca = PCA(pcs)
#projected = pca.fit_transform(np.log10(enrichments_pca.values.transpose()+1))
projected = pca.fit_transform(enrichments_pca.values.transpose())

#s_table = copy.deepcopy(s_table)
s_table = s_table[s_table["sample_group"] != "healthy control"]
#p_table = copy.deepcopy(p_table)

pc_scaled = []
pc_proj = []
for pc in range(pcs):
    pc_proj.append(f"PC-{pc+1}-proj")
    s_table[f"PC-{pc+1}-proj"] = projected[:, pc]
    
    ys = projected[:, pc]
    scaley = 1.0/(ys.max() - ys.min())
    s_table[f"PC-{pc+1}-scaled"] = ys * scaley
    pc_scaled.append(f"PC-{pc+1}-scaled")

for pc in range(pcs):
    p_table[f"loading-{pc+1}"] = pca.components_[pc, :]


largest_proj = max(s_table["PC-3-proj"])

fig = plt.figure(constrained_layout=True, figsize=[12, 10])
axd = fig.subplot_mosaic(
    """
    CCDDH
    CCDDH
    AABBI
    AABBI
    JJJJJ
    EEEEE
    FFFFF
    GGGGG
    """,
    gridspec_kw={
        "wspace": 0.75,
        "hspace": 0.85,
    },
)

for epitope, metadata in EPITOPES.items():
    if epitope in ["CTD-1", "CTD-2", "CTD-3"]: continue
    limits = metadata["limits"]
    pep_id = p_table.loc[p_table["Loc"].isin(range(limits[0], limits[1])),:].index
    p_table.loc[pep_id, "epitope"] = epitope

s_table["sample_color"] = [
        sample_group_colors[sg] 
        for sg in s_table["sample_group"]
]

for group, group_s in s_table.groupby("sample_group"):
    for subp, p1, p2 in zip(["C", "D"], (1, 1), (2, 3)):
        axd[subp].scatter(
            x = group_s[f"PC-{p1}-scaled"], 
            y = group_s[f"PC-{p2}-scaled"],
            c = group_s["sample_color"],
            label = group,
            s=35,
            alpha=0.8,
            edgecolors="black"
        )
        explained_variance_ratio_ = np.zeros(10)

for row, peptide in p_table.iterrows():
    c = epitope_colors[peptide["epitope"]]
    a = 0.8 if peptide["epitope"] != "none" else 0.3
    for subp, p1, p2 in zip(["A", "B"], (1, 1), (2, 3)):
        axd[subp].arrow(
            0,0, peptide[f"loading-{p1}"], peptide[f"loading-{p2}"], 
            color=c, width=0.002, alpha=a, edgecolor="black"
        )

p_table.to_csv("./pca_p_table.csv")

for i, lo in enumerate(["E", "F", "G"]):
    for epitope, metadata in EPITOPES.items():
        if epitope in ["CTD-1", "CTD-2", "CTD-3"]: continue
        e_l = metadata["limits"]
        axd[lo].plot(
            p_table["Loc"].astype(int).values[e_l[0]:e_l[1]],
            p_table[f"loading-{i+1}"].values[e_l[0]:e_l[1]],
            c = epitope_colors[epitope],
            lw = 4,
            label=f"{e} - [{e_l[0]}, {e_l[1]})"
        )

    axd[lo].plot(
        p_table["Loc"].astype(int),
        p_table[f"loading-{i+1}"],
        c = "black",
        lw = 0.7
    )

axd["J"].plot(
    p_table["Loc"].astype(int),
    np.repeat(-10, len(p_table)),
    c = "black",
    lw = 0.7
)
axd["J"].set_ylim(-5, 5)
axd["J"].axis('off')

protein_regions = {
    "NTD" : [29, 294],
    "RBD" : [332, 523],
    "CTD" : [530, 681],
    "FP" : [816, 833],
    "HR1" : [910, 987],
    "HR2" : [1162, 1204]
}

for pr, limits in protein_regions.items():

    start = limits[0]
    width = (limits[1] - limits[0]) 
    middle = start + (width//2)
    anno = f"{pr}\n[{limits[0]},\n {limits[1]})"

    rect_v = patches.Rectangle(
            (start, -10),
            width=width,
            height=3,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor="white"
    )

    axd["J"].text(middle, -2, anno, va="center", ha="center", size=12)
    axd["J"].add_patch(rect_v)

for ep, metadata in EPITOPES.items():
    if ep in ["CTD-1", "CTD-2", "CTD-3"]: continue
    limits = metadata["limits"]

    start = limits[0]
    width = (limits[1] - limits[0]) 
    middle = start + (width//2)

    rect_v = patches.Rectangle(
            (start, -10),
            width=width,
            height=1,
            clip_on=False, 
            linewidth=1, 
            edgecolor='black',
            facecolor=epitope_colors[ep]
    )

    axd["J"].add_patch(rect_v)

epitope_order = [
    "NTD",
    "CTD",
    "FP", 
    "HR2",
]

axd["I"].axis('off')
legend_elements = [
    patches.Patch(
         facecolor=epitope_colors[epit], 
         edgecolor="black", 
         label=epit
    )
    for i, epit in enumerate(epitope_order)
]

legend_elements.append(
    patches.Patch(
         facecolor="black", 
         edgecolor="black", 
         label="None"
    )
)

axd["I"].legend(
        handles=legend_elements, 
        loc="lower left",
        bbox_to_anchor = [-0.5, 0.0],
        frameon=False
)

axd["H"].axis('off')
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
        loc="lower left",
        bbox_to_anchor = [-0.5, 0.0],
        frameon=False
)


kw = dict(ha="center", va="center", fontsize=15, color="black")
axd["C"].text(-0.05, 1.10, f"A", transform=axd["C"].transAxes, **kw)
axd["A"].text(-0.05, 1.10, f"B", transform=axd["A"].transAxes, **kw)
axd["E"].text(-0.05, 1.10, f"C", transform=axd["E"].transAxes, **kw)

axd["C"].set_xlabel(
    f"PC 1 ({round(pca.explained_variance_ratio_[0]*100, 2)}% EVR)"
)
axd["C"].set_ylabel(
    f"PC 2 ({round(pca.explained_variance_ratio_[1]*100, 2)}% EVR)"
)

axd["D"].set_xlabel(
    f"PC 1 ({round(pca.explained_variance_ratio_[0]*100, 2)}% EVR)"
)
axd["D"].set_ylabel(
    f"PC 3 ({round(pca.explained_variance_ratio_[2]*100, 2)}% EVR)"
)

axd["A"].set_xlabel("PC 1")
axd["A"].set_ylabel("Loadings\nPC 2")
axd["C"].set_ylabel("Scores\nPC 2")
axd["B"].set_xlabel("PC 1")
axd["B"].set_ylabel("PC 3")

plt.tight_layout()
fig.savefig(sys.argv[2])

