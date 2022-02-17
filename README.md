# Comprehensive characterization of the antibody responses to SARS-CoV-2 Spike protein finds additional vaccine- induced epitopes beyond those for mild infection
**Published in** [elife](http://dx.doi.org/10.7554/eLife.73490).
Authors: Meghan E. Garrett\*, Jared G. Galloway\*, Caitlin Wolf, Jennifer K. Logue, Nicholas Franko, Helen Chu, Frederick A. Matsen IV^, Julie Overbaugh^

\* these authors contributed equally to this work.
^ co-corresponding authors


### What is this?

This repository is here to serve as a static archive for all analysis done in 
*Comprehensive characterization of the 
antibody responses to SARS- CoV- 2 Spike 
protein finds additional vaccine- induced 
epitopes beyond those for mild infection* -- an elife published manuscript which
may be found at 
[http://dx.doi.org/10.7554/eLife.73490](http://dx.doi.org/10.7554/eLife.73490).

Primarily, this repository aims to provide transparency and reproducibility within the context of our study. 
The files provide a complete set of materials to replicate the analysis (from fastq to figures) 
found in our manuscript with a single execution of a Nextflow pipeline.
For exploring the data interactively, please see our
[DMS-View data repository](https://github.com/matsengrp/vacc-dms-view-host-repo). 

Ultimately, running the pipeline will result in an 
[xarray DataSet](http://xarray.pydata.org/en/stable/), 
(see [phippery](https://github.com/matsengrp/phippery), 
for more on this dataset organization), 
as well as Figures as seen in the manuscript.
The pipeline runs the analysis and plotting code for two sets of of phage-display library batch replicates. The figure sets for each respective batch are separated here:

1. Figure set from ["SPIKE1 Replicates"](Manuscript-Figures/SPIKE1/) 
2. Figure set from ["SPIKE2 Replicates"](Manuscript-Figured/SPIKE2/) <- presented in the manuscript

If interested in obtaining raw data to perform the analysis yourself,
please feel free to contact Jared Galloway:
jgallowa (at) fredhutch (dot) org. 


### Exploring the data, interactively

Running the pipeline here is **not** the suggested approach for exploring our data.
While running the pipeline is quite simple with some configuration (see *Running the Pipeline*),
it involves processing over 600 sequence alignments and running downstream
esoteric analysis/plotting code specific to our sample's 
[metadata](nextflow-pipeline-config/sample_table.csv). Thus, tweaking parameters may be a headache.

Instead, if you're interested in simply exploring our the rich amount of data from this study,
we strongly suggest checking out the pre-processed and publicly explorable 
[DMS-View data repository](https://github.com/matsengrp/vacc-dms-view-host-repo). 
There, we have formatted and hosted the data for every sample in the study
(398 replicates across two library batches of phage display) 
to be used with the amazing DMS-View tool put out by the
[Bloom Lab](https://research.fredhutch.org/bloom/en.html). 
For more on this, see the repository
[README](https://github.com/matsengrp/vacc-dms-view-host-repo/blob/main/README.md)

### Material Overview

We provide a fully reproducible and automated workflow which ingests 
raw sequencing data and performs all analyses presented in the paper. 
The workflow extends our more generalized PhIP-Seq alignment pipeline, 
[PhIP-flow](https://github.com/matsengrp/phip-flow)
The materials for analysis are primarily broken down into three categories:

1. `image-template/` The configuration scripts defining a container image, which is used to build 
        the container with all version-specific [phippery](https://github.com/matsengrp/phippery) source code along with other non-local python package dependencies for analysis and plotting.
        
2. `analysis-scripts/` The python scripts for computing normalizations on the data, as well as plotting code to produce our final figures. 

3. `nextflow-pipeline-config/` The Nextflow pipeline script as well as all necessary configuration scripts to run the workflow either (a) locally on a computer with docker installed, or (b) a [SLURM](https://slurm.schedmd.com/documentation.html) managed cluster with singularity available. 


### Running the Pipeline

*What do I need?* 

Working installation of 
[Docker](https://docs.docker.com/get-docker/) and 
[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). 
Maybe some computing power if starting from raw fastq.

*How do I run it?* 

1.  For running locally (not recommended) install Docker + Nextflow. Otherwise, 
we have a configuration script that would take very little editing to run the analysis on a [SLURM](https://slurm.schedmd.com/documentation.html) managed cluster with access to Nextflow and Singularity modules

2. Clone this repository and obtain the raw fastq sequences. The data will likely come in the form of a tarball archive which when extracted, will provide an `NGS/` folder containing all the demultiplexed sample sequence data as described in our [sample table](nextflow-pipeline-config/sample_table.csv). Place the `NGS/` directory within the repository's [nextflow-pipeline-config](./nextflow-pipeline-config) subdirectory.

3. Generate a config script specific to your compute infrastructure. Consult the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on fitting the parameters to your specific infrastructure. 
We provide an example of such a configuration for our Fred Hutch SLURM managed cluster in 
[this file](./nextflow-pipeline-config/phipflow.config.bt2). 
Whatever your configuration, the file *must* include parameters specified in the `PARAMS{...}` block of our config script.

4. Run the pipeline. An example of how we call the `nextflow run` command on our compute infrastructure the pipeline can be seen in the `nextflow-pipeline-config/run_analysis.sh`

```bash
(base) quokka phage-dms-vacc-analysis/nextflow-pipeline-config ‹master*› » ./run_analysis.sh 
N E X T F L O W  ~  version 20.07.1
Launching `PhIP-analysis.nf` [golden_ekeblad] - revision: 02870c3fbe
[01/807cac] process > generate_fasta_reference (1) [100%] 1 of 1, cached: 1 ✔
[f3/915808] process > generate_index (1)           [100%] 1 of 1, cached: 1 ✔
[95/1f1df0] process > short_read_alignment (633)   [100%] 633 of 633, cached: 633 ✔
[25/3262fa] process > sam_to_stats (633)           [100%] 633 of 633, cached: 633 ✔
[40/d53f0f] process > sam_to_counts (633)          [100%] 633 of 633, cached: 633 ✔
[37/3ab3b4] process > collect_phip_data (1)        [100%] 1 of 1, cached: 1 ✔
[22/efe9a0] process > compute_enrichment_stats (1) [100%] 1 of 1, cached: 1 ✔
[2d/f9fe6f] process > analysis_plotting (1)        [100%] 1 of 1, cached: 1 ✔

71.45user 5.57system 0:20.03elapsed 384%CPU (0avgtext+0avgdata 1743320maxresident)k
8200inputs+39768outputs (4major+906688minor)pagefaults 0swaps
```

The pipeline will put all batch-specific figures and the respective xarray datasets in the `phip_data_dir` as defined by the phip-flow configuration scripts. 

### Static containers

vacc-ms-analysis:vacc-ms-analysis [![vacc-ms-analysis](https://quay.io/repository/matsengrp/vacc-ms-analysis/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/vacc-ms-analysis) An extension of the [phippery](https://github.com/matsengrp/phippery) container with all back end function source code and dependencies listed in the `image-template/requirements.txt`

phippery:vacc-ms-analysis [![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery) phippery container

quay.io/jgallowa/bowtie2:vacc-ms-analysis - 
a static container containing the
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment tool. 
The original image was hosted by [Biocontainers](https://biocontainers.pro/)

quay.io/matsengrp/samtools-1.3:vacc-ms-analysis
a static container containing the
[samtools](http://www.htslib.org/) software. 
The original image was hosted by [Biocontainers](https://biocontainers.pro/)
