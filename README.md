# Differences in Spike epitope targeting between antibodies elicited by SARS-CoV-2 mRNA vaccines versus infection. 

### What is this?

This repository contains a complete set of materials, excluding raw data, 
to replicate the analysis found in our manuscript.
If interested in obtaining raw data,
please feel free to contact [jared galloway]():
jgallowa (at) fredhutch (dot) org. 

## Abstract

Control of the COVID-19 pandemic will rely on SARS-CoV-2 vaccine-elicited antibodies to protect against emerging and future variants; a more complete description of the differential humoral responses from infection or vaccination is needed. Here we comprehensively profiled the linear epitopes and pathways of escape for Spike-specific antibodies from individuals either enrolled in a phase 1 clinical trial of the mRNA-1273 Moderna vaccine (n=49) or enrolled in a Seattle-based study of people with diverse SARS-CoV-2 infection and/or vaccination status (n=60). Four epitopes accounted for most of the variance between samples: N-terminal domain (NTD), C-terminal domain (CTD), fusion peptide (FP), and heptad repeat 2 (HR2) epitopes. Binding to the FP and HR2 epitopes alone was associated with mild infection, whereas those with severe infection or vaccination had antibodies that bound to all four epitopes. Epitope binding appeared to change over time after vaccination, but other covariates such as mRNA vaccine dose, vaccine type (Pfizer BNT162b2 or Moderna mRNA-1273), or participant age did not appear to affect antibody binding to these epitopes. Vaccination induced a strong, uniform escape profile in the NTD, CTD, and HR2 regions, whereas infection elicited a strong response in the FP region with an antibody escape profile that was maintained after vaccination. Overall, vaccination led to a greater number of epitopes targeted across Spike and induced a uniform escape profile across individuals in many cases. These findings have implications for the selection of SARS-CoV-2 escape mutations on a population level. 


### overview


We provide a fully reproducible automated workflow which ingests raw sequencing data and performs all analyses presented in the paper. 
The workflow defines and runs the processing steps within publicly available and static Docker software containers, 
including phippery and phip-flow described in the Methods section. 
The source code, Nextflow script, software dependencies, and instructions for re-running the analysis 
The materials for analysis are primarily broken down into three categories:

1. `image-template/` The configuration scripts defining a container image, which is used to build 
        a container with all version-specific [phippery](https://github.com/matsengrp/phippery) source code along with other non-local python package dependencies.
        
vacc-ms-analysis: [![vacc-ms-analysis](https://quay.io/repository/matsengrp/vacc-ms-analysis/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/vacc-ms-analysis)

2. `analysis scripts/` the python scripts which perform analysis given a set of parameters followed by generation of all parameter-specific plots as those seen in the manuscript.  




### What do I need to run this?

Docker and Nextflow. Maybe some computing power if starting from raw fastq.

1. install docker and nextflow command line tools, I'm using the following versions

```
(phippery) ubuntu analysis/main-figs-scripts ‹master*› » docker -v
Docker version 20.10.1, build 831ebea
(phippery) ubuntu analysis/main-figs-scripts ‹master*› » nextflow -v
nextflow version 20.04.1.5335
```

## Figure PDF's

## Sample, Peptide and Epitope Binding Region Metadata.

## PhIP-Flow Sample Alignment Pipeline

## Phippery Downstream Analysis Scripts and Dependencies

The downstream analysis described in the manuscript relies on an in-house open-source python 
package, [Phippery](), as well as some other popular python packages. 
These version-specific dependencies are all contained in the container images 
listed in the `pipeline-scripts/phipflow.config` 

## Can I run this on my own data?

short answer, yes, but the configuration, downstream analysis, and plotting
code are all taylored to the `pipeline-scripts/sample_table.csv`. 
At the time of publication, interactive visualization tools for many of the analysis
seen here are currently being worked on at 

The steps involved in using this for your own data would involve creating your own sample metadata. We're in the process of creating interactive apps for many of the analysis and vidsualizations here. 

So ... here's basically my best description. Again, feel free to reach out with any questions

## 


