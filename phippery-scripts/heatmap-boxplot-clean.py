
import phippery
#from phippery.normalize import svd_rank_reduction,enrichment, standardized_enrichment, differential_selection_wt_mut, rank_data, svd_aa_loc, size_factors, cpm
from phippery.collapse import collapse_sample_groups, pairwise_correlation_by_sample_group
from phippery.utils import *
from phippery.tidy import tidy_ds
from phippery.modeling import neg_binom_model
import xarray as xr
import numpy as np
import pandas as pd
import pickle
from collections import defaultdict
import copy
from plotnine import *
from numpy.linalg import svd
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import to_hex
import matplotlib.patches as patches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
import seaborn as sns
import time
import warnings
from adjustText import adjust_text
warnings.filterwarnings("ignore")
import numpy as np
from sklearn.decomposition import PCA
from statannot import add_stat_annotation
import sys

######
#KNOBS
######

#metric = "standardized_enrichment"
metric = "counts_enrichment"
#metric = "neg_binom_mlxp"
highlight="sample_group"
batch = "SPIKE1"
np.random.seed(24)
mode_agg="wt_only"
outdir="heatmap-boxplot/"
ag="sum"
sample_colormap="Set3"
#epitope_colormap="Accent"
epitope_colormap="Set1"
heatcolormap="RdPu"
window_size=10
pcs=3


ds = phippery.load(sys.argv[1])
if batch in ["SPIKE1", "SPIKE2"]:
    batch_samples = id_coordinate_subset(ds, where="library_batch", is_equal_to=batch)
    ds = ds.loc[dict(sample_id=batch_samples)]

if not os.path.exists(outdir): os.mkdir(outdir)

# Helper function used for visualization in the following examples
def identify_axes(ax_dict, fontsize=12):
    """
    Helper to identify the Axes in the examples below.
    Draws the label in a large font in the center of the Axes.
    Parameters
    ----------
    ax_dict : dict[str, Axes]
        Mapping between the title / label and the Axes.
    fontsize : int, optional
        How big the label should be.
    """
    kw = dict(ha="center", va="center", fontsize=fontsize, color="black")
    for k, ax in ax_dict.items():
        if k not in ["A", "B"]: continue
        ax.text(-0.05, 1.10, f"{k})", transform=ax.transAxes, **kw)

# transform colormap
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


hcmap = shiftedColorMap(getattr(matplotlib.cm, heatcolormap), midpoint=0.4)
#hcmap = getattr(matplotlib.cm, heatcolormap)

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
    "HR2" : [1135,1165],
}

cmap=plt.get_cmap(epitope_colormap)
epitope_colors = {
    "none":"black",
}
for i, (e, _) in enumerate(epitope_limits.items()):
    epitope_colors[e] = cmap.colors[7]

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

print(batch_beads)

wt = id_coordinate_subset(ds, where="is_wt", is_equal_to=True, table="peptide_table")
ds = ds.loc[dict(peptide_id=wt)]
enr_label = "Enrichment"

"""
if metric == "counts_enrichment":

    ds = ds.loc[dict(peptide_id=wt, sample_id=list(set.union(batch_lib, batch_emp)))]
    enrichment(ds, list(batch_lib), "counts", new_table_name="counts_enrichment")

    ds = ds.loc[dict(sample_id=list(batch_emp))]
    enr_label = "Enrichment"

elif metric == "standardized_enrichment":

    ds = ds.loc[dict(peptide_id=wt, sample_id=list(set.union(batch_lib, batch_emp, batch_beads)))]
    standardized_enrichment(
            ds, 
            list(batch_lib), 
            list(batch_beads), 
            "counts", 
            new_table_name="standardized_enrichment"
    )

    ds = ds.loc[dict(sample_id=list(batch_emp))]
    enr_label = "Standardized Fold enrichment"

else:

    ds = ds.loc[dict(peptide_id=wt, sample_id=list(set.union(batch_beads, batch_emp)))]
    size_factors(ds)

    beads_ds = ds.loc[dict(sample_id=list(batch_beads))]
    ds = ds.loc[dict(sample_id=list(batch_emp))]

    neg_binom_model(
        ds,
        beads_ds,
        nb_p=2,
        trim_percentile=99.9,
        data_table="size_factors",
        inplace=True,
        new_table_name="neg_binom_mlxp"
    )
    enr_label = "Negative Binomial MLXP"
"""
############################################
################ CURATE SAMPLES ############
############################################

# remove 250 dose samples
#mg_250 = set(id_coordinate_subset(
#    ds, 
#    where="dose", 
#    is_equal_to="250ug"
#))

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

#cohort_250 = set.intersection(
#        mg_250,
#        cohort
#    )
cohort_119 = set.intersection(
        time_point,
        cohort
    )

print(len(cohort_119))

sams = set(ds.sample_id.values)
#sams = sams - set.union(cohort_119, cohort_250)
sams = sams - cohort_119
#sams = sams - set.unioncohort_250)
ds = ds.loc[dict(sample_id=list(sams))]

############################################
################ PREP DATA #################
############################################

# Grab enrichments
enrichments = ds[f"{metric}"]

#sys.exit()

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


#identify_axes(axd)

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

boxpairs_dict = {
    "B" : [('Vaccinated - (No prior infection) Cohort 1', 'Non-Hospitalized Serum Cohort 2')],
    "C" : [('Vaccinated - (No prior infection) Cohort 1', 'Non-Hospitalized Serum Cohort 2')],
    "D" : [('Vaccinated - (No prior infection) Cohort 1', 'Non-Hospitalized Serum Cohort 2')],
    "E" : [('Vaccinated - (No prior infection) Cohort 1', 'Non-Hospitalized Serum Cohort 2')],
}

boxpairs = [
        ("Non-Hospitalized Serum Cohort 2",g) 
        for g in set(sample_group_order)-set(["Non-Hospitalized Serum Cohort 2"])
]

first=True
for ax, epitope in zip(["B", "C", "D", "E"], epitope_order):
    sns.boxplot(
            data=bxr_s_table, 
            x="sample_group", 
            y=f"{epitope}_raw_agg", 
            ax=axd[ax], 
            palette=sample_group_colors,
    )

    add_stat_annotation(
            axd[ax], 
            data=bxr_s_table, 
            x="sample_group", 
            y=f"{epitope}_raw_agg", 
            order=sample_group_order,
            #box_pairs = boxpairs_dict[ax],
            box_pairs = boxpairs,
            test='Mann-Whitney', 
            comparisons_correction=None,
            text_format='star', 
            loc='inside', 
            verbose=2
    )
    limits = epitope_limits[epitope]

    ylabel = "Summed enrichment"
    title = f"{epitope} Epitope \n{limits[0]} - {limits[1]}"
    #ylabel = "$\sum_{WT_{i}}$" + f" s.t. $i \in [{limits[0]},{limits[1]})$"
    axd[ax].set_ylabel(ylabel)
    axd[ax].set_xlabel("")
    hue_color = [min(1, i*0.6) for i in epitope_colors[epitope]]
    axd[ax].set_title(title, fontdict={'color':'slategray'})
    axd[ax].get_xaxis().set_visible(False)
    first = False


############################################
################ HEATMAP ###################
############################################

hm_s_table = copy.deepcopy(s_table)
hm_p_table = copy.deepcopy(p_table)

windows = list(range(0, max(hm_p_table["Loc"]), window_size))
enrichment_columns = []
for l in range(len(windows)-1):
    start = windows[l]
    end = windows[l+1]
    pep_id = hm_p_table.loc[hm_p_table["Loc"].isin(range(start, end)),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]
    enrichment_columns.append(f"[{start}, {end})")
    if ag == "sum":
        agg = epitope_enrichment.sum(axis=0).values
    else:
        agg = epitope_enrichment.mean(axis=0).values
    hm_s_table[f"[{start}, {end})"] = [min(150, v) for v in agg]


sample_order = []
for g, g_df in hm_s_table.groupby(
        ["vaccination_status", "manuscript_cohort", "sample_group"]
        ):
    sample_order.extend(g_df.index)
    
hm_s_table = hm_s_table.reindex(sample_order)

s_e = hm_s_table[enrichment_columns]

cbar_kws = dict(use_gridspec=False, location="right", label=enr_label)
sns.heatmap(s_e, ax = axd["A"], cmap=hcmap, cbar_kws=cbar_kws)
axd["A"].get_yaxis().set_visible(False)
axd["A"].set_ylabel("Amino acid position windows")


base=0
for g, g_df in hm_s_table.groupby(["manuscript_cohort", "sample_group"]):
    
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
for g, g_df in hm_s_table.groupby(["manuscript_cohort"]):
    
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
for g, g_df in hm_s_table.groupby(["vaccination_status"]):
    
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

for ep, limits in epitope_limits.items():

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
            #facecolor=epitope_colors[ep]
            facecolor="slategrey"
    )

    axd["A"].add_patch(rect_v)

axd["A"].get_xaxis().set_visible(False)

legend_elements = [
    Patch(
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
        #loc="upper right"
         ha="right", va="top",
         #bbox=dict(boxstyle="square",
         #          ec=(1., 0.5, 0.5),
         #          fc=(1., 0.8, 0.8),
         #          )
         )

"""
pca_p_table = pd.read_csv("../pre-processed-datasets/pca_p_table.csv")
for i, lo in enumerate(["F", "I", "J"]):
    for e, e_l in epitope_limits.items():
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
"""

axd["F"].get_xaxis().set_visible(False)
axd["F"].yaxis.set_label_position("right")
axd["F"].set_ylabel("Comp 1\nLoadings")

axd["I"].get_xaxis().set_visible(False)
axd["I"].yaxis.set_label_position("right")
axd["I"].set_ylabel("Comp 2\nLoadings")

axd["J"].set_xlabel("Amino Acid Position")
axd["J"].yaxis.set_label_position("right")
axd["J"].set_ylabel("Comp 3\nLoadings")
axd["J"].yaxis.label.set_fontsize(10)

axd["C"].set_ylabel("")
axd["E"].set_ylabel("")


for a in axd:
    ax = axd[a]
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(14)

"""
legend_elements = [
    Patch(
         facecolor=epitope_colors[ep], 
         edgecolor="black", 
         label=ep + " Epitope region"
    )
    for ep in epitope_order
]

axd["H"].legend(
        handles=legend_elements, 
        loc="upper right",
        frameon=False
)
"""
plt.tight_layout()
#fig.savefig(f"{outdir}/{metric}-{batch}-36d-grey.png")
fig.savefig(sys.argv[2])
