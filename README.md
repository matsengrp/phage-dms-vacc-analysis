# Differences in Spike epitope targeting between antibodies elicited by SARS-CoV-2 mRNA vaccines versus infection. 

### What is this?

This repository contains a complete set of materials, excluding raw data, 
to replicate the analysis found in our manuscript.
If interested in obtaining raw data,
please feel free to contact [jared galloway]():
jgallowa (at) fredhutch (dot) org. 

## Abstract

Control of the COVID-19 pandemic will rely on SARS-CoV-2 vaccine-elicited antibodies to protect against emerging and future variants; a more complete description of the differential humoral responses from infection or vaccination is needed. Here we comprehensively profiled the linear epitopes and pathways of escape for Spike-specific antibodies from individuals either enrolled in a phase 1 clinical trial of the mRNA-1273 Moderna vaccine (n=49) or enrolled in a Seattle-based study of people with diverse SARS-CoV-2 infection and/or vaccination status (n=60). Four epitopes accounted for most of the variance between samples: N-terminal domain (NTD), C-terminal domain (CTD), fusion peptide (FP), and heptad repeat 2 (HR2) epitopes. Binding to the FP and HR2 epitopes alone was associated with mild infection, whereas those with severe infection or vaccination had antibodies that bound to all four epitopes. Epitope binding appeared to change over time after vaccination, but other covariates such as mRNA vaccine dose, vaccine type (Pfizer BNT162b2 or Moderna mRNA-1273), or participant age did not appear to affect antibody binding to these epitopes. Vaccination induced a strong, uniform escape profile in the NTD, CTD, and HR2 regions, whereas infection elicited a strong response in the FP region with an antibody escape profile that was maintained after vaccination. Overall, vaccination led to a greater number of epitopes targeted across Spike and induced a uniform escape profile across individuals in many cases. These findings have implications for the selection of SARS-CoV-2 escape mutations on a population level. 


### Overview

We provide a fully reproducible automated workflow which ingests raw sequencing data and performs all analyses presented in the paper. 
The workflow defines and runs the processing steps within publicly available and static Docker software containers, 
including phippery and phip-flow described in the Methods section. 
The source code, Nextflow script, software dependencies, and instructions for re-running the analysis 
The materials for analysis are primarily broken down into three categories:

1. `image-template/` The configuration scripts defining a container image, which is used to build 
        a container with all version-specific [phippery](https://github.com/matsengrp/phippery) source code along with other non-local python package dependencies.
        
2. `analysis-scripts/` the python scripts which perform analysis given a set of parameters followed by generation of all parameter-specific plots as those seen in the manuscript.  

3. `nextflow-pipeline-config` All the necessary configuration scripts to run the pipeline either locally on a computer with docker installed, or a [SLURM](https://slurm.schedmd.com/documentation.html) managed cluster with singularity available. 



### Quickstart

Docker and Nextflow. Maybe some computing power if starting from raw fastq.

1. For running locally (not reccomended) install Docker + Nextflow. Otherwise,
we have a configuration script that would take very little editing to run on a [SLURM](https://slurm.schedmd.com/documentation.html managed cluser with access to Nextflow and Singularity modules)

2. Clone this repository and obtain the raw fastq sequences -- being sure to put them in the nextflow-pipeline-config directory under in the subdirectory names `NGS/`. 

3. Inside the `nextflow-pipeline-config/phipflow.config.bt2` script, you may setup configuration settings for your particular computing environment i.e. which partitions get used for running each of the jobs as well as the resources allocated. For more information about setting up the configuration for your machine, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html).

4. Run the pipeline!

```bash

```



### Output and analysis parameters



### Static containers

vacc-ms-analysis:vacc-ms-analysis [![vacc-ms-analysis](https://quay.io/repository/matsengrp/vacc-ms-analysis/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/vacc-ms-analysis) An extension of the phippery container with all backend function source code and dependencies listed in the `image-template/requirements.txt`

phippery:vacc-ms-analysis [![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery) phippery container

quay.io/jgallowa/bowtie2:vacc-ms-analysis

quay.io/matsengrp/samtools-1.3:vacc-ms-analysis
