/*
PhIP-Flow Pipeline Script

Author: Jared G. Galloway
*/


peptide_reference_ch = Channel.fromPath("${params.peptide_table}")

// CONVERT PEPTIDE METADATA TO FASTA
process generate_fasta_reference {

    label 'single_thread_large_mem'

    input: file "pep_ref" from peptide_reference_ch

    output: file "peptides.fasta" into pep_channel_fasta

    shell:    
        """
        phippery peptide-md-to-fasta -d ${pep_ref} -o peptides.fasta
        """
}


// GENERATE INDEX
process generate_index {
 
    label 'multithread'

    input:
        file "pep_fasta" from pep_channel_fasta

    output:
        set(
            val("peptide_ref"), 
            file("peptide_index") 
        ) into pep_channel_index

    shell:    

        if ("$params.alignment_tool" == "bowtie")
            """
            mkdir peptide_index
            bowtie-build --threads 4 \
            ${pep_fasta} peptide_index/peptide
            """

        else if ("$params.alignment_tool" == "bowtie2")

            """
            mkdir peptide_index
            bowtie2-build --threads 4 \
            ${pep_fasta} peptide_index/peptide
            """
}


// CREATE SAMPLE CHANNEL
Channel
    .fromPath("${params.sample_table}")
    .splitCsv(header:true)
    .map{ row -> 
        tuple(
            "peptide_ref",
            row.sample_id,
            new File("${row.seq_dir}/${row.fastq_filename}")
        ) 
    }
    .set { samples_ch }

index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            sample[1],          // sample_id
            file(ref[1]),       // index files
            file(sample[2]),    // sample path
        )
    }


// ALIGN ALL SAMPLES TO THE REFERENCE
process short_read_alignment {

    label 'multithread'
    //container = 'quay.io/jgallowa/bowtie2:vacc-ms-analysis' 

    input:
        set( 
            val(sample_id), 
            file(index),
            file(respective_replicate_path),
        ) from index_sample_ch

    output:
        set(
            val(sample_id),
            file("${sample_id}.sam")
        ) into aligned_reads_sam

    shell:

        if ("$params.alignment_tool" == "bowtie")
            """
            ${params.fastq_stream_func} ${respective_replicate_path} | \
            bowtie ${params.align_args} --sam -x ${index}/peptide - > ${sample_id}.sam
            """

        else if ("$params.alignment_tool" == "bowtie2")
            """
            ${params.fastq_stream_func} ${respective_replicate_path} | \
            bowtie2 ${params.align_args} -x ${index}/peptide - > ${sample_id}.sam
            """
}


// SPLIT CHANNEL FOR COUNTS AND STATS IN PARALLEL
aligned_reads_sam.into{aligned_reads_for_counts; aligned_reads_for_stats}


// COMPUTE ALIGNMENT STATS FOR ALL STATS
process sam_to_stats {

    label 'multithread'
    ////container = 'quay.io/biocontainers/samtools:1.3--h0592bc0_3'
    //container = 'quay.io/matsengrp/samtools-1.3:vacc-ms-analysis'

    input:
        set(
            val(sample_id),
            file(sam_file)
        ) from aligned_reads_for_stats
    
    output:
        file("${sample_id}.stats") into alignment_stats_ch
    
    shell:
        """
        samtools stats ${sam_file} | grep ^SN | cut -f 2- | 
        sed '1p;7p;22p;25p;d' > ${sample_id}.stats
        """ 
}


// COMPUTE COUNTS FOR ALL SAMPLES
process sam_to_counts {
    
    label 'multithread'
    //container = 'quay.io/matsengrp/samtools-1.3:vacc-ms-analysis'

    input:
        set(
            val(sample_id),
            file(sam_file)
        ) from aligned_reads_for_counts

    output:
        file("${sample_id}.counts") into counts_ch

    script:
        """
        samtools view -u -@ 28 ${sam_file} | \
        samtools sort -@ 28 - > ${sample_id}.bam
        samtools sort -@ 28 ${sample_id}.bam -o ${sample_id}.sorted 
        mv ${sample_id}.sorted ${sample_id}.bam
        samtools index -b ${sample_id}.bam
        samtools idxstats ${sample_id}.bam | \
        cut -f 1,3 | sed "/^*/d" > ${sample_id}.counts
        """
}

// COLLECT AND MERGE ALL 
process collect_phip_data {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}", mode: 'copy'
    label 'single_thread_large_mem'
    //container = 'quay.io/matsengrp/phippery:vacc-ms-analysis' 

    input:
        file all_counts_files from counts_ch.collect()
        file all_alignment_stats from alignment_stats_ch.collect()
        file sample_table from Channel.fromPath("${params.sample_table}")
        file peptide_table from Channel.fromPath("${params.peptide_table}")

    output:
        file "${params.dataset_prefix}.phip" into phip_data_ch

    script:
        """
        phippery collect-phip-data \
        -s_meta ${sample_table} \
        -p_meta ${peptide_table} \
        -c '*.counts' \
        -s '*.stats' \
        -o ${params.dataset_prefix}.phip
        """ 
}

// RUN ALL ANALYSIS
process compute_enrichment_stats {
    
    //publishDir "${params.phip_data_dir}/", mode: 'copy'
    publishDir "${params.phip_data_dir}/${params.alignment_tool}", mode: 'copy'
    label 'single_thread_large_mem'
    ////container = 'quay.io/matsengrp/phippery:latest' 
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 

    input:
        file phip_ds from phip_data_ch
        file alignment_stats from Channel.fromPath("../analysis-scripts/alignment-stats.py")
        file analysis from Channel.fromPath("../analysis-scripts/layer-enrichment-stats.py")
        

    output:
        file "alignment-stats.pdf"
        file "layered-analysis.phip" into layered_phip_data_ch

    """
    set -eu
    python  ${alignment_stats} -dataset ${phip_ds} -out alignment-stats.pdf
    python ${analysis} -dataset ${phip_ds} -out layered-analysis.phip
    """ 
}

batches = Channel.from("SPIKE2", "SPIKE1")
batch_layer = batches.combine(layered_phip_data_ch)
epitopes = Channel.fromPath("../analysis_scripts/epitopes.py")
bde = batch_layer.combine(epitopes).into(6)

pca = Channel.fromPath("../analysis-scripts/pca-scatter-directions.py").combine(bde[0])
hmb = Channel.fromPath("../analysis-scripts/heatmap-boxplot.py").combine(pca)

// RUN PCA and Plot heatmap
process pca_heatmap {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(heatmap),
            file(pca),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from hmb

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${pca} -dataset ${layered_phip_ds} -batch ${batch} -out pca.${params.plot_format}
        python ${heatmap} -dataset ${layered_phip_ds} -batch ${batch} -out  heatmap-boxplot.${params.plot_format}
        """ 
}

logo = Channel.fromPath("../analysis-scripts/logopairs.py").combine(bde[1]).into(2)

// RUN logoplot 1
process logoplots_moderna {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(logo),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from logo[0]

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${logo} -dataset ${layered_phip_ds} -subgroup moderna -batch ${batch} -out logopairs-moderna.${params.plot_format}
        """ 
}

// RUN logoplot 2
process logoplots_haarvi {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(logo),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from logo[1]

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${logo} -dataset ${layered_phip_ds} -subgroup haarvi -batch ${batch} -out logopairs-haarvi.${params.plot_format}
        """ 
}

thr = Channel.fromPath("../analysis-scripts/threshold-epi-binding-barplot.py").combine(bde[2]).into(2)

process threshold_haarvi {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(thresh),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from thr[0]

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${thresh} -dataset ${layered_phip_ds} -subgroup haarvi -out epitope_wt_thresholds-haarvi.${params.plot_format}
        """ 
}

process threshold_moderna {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(thresh),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from thr[1]

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${thresh} -dataset ${layered_phip_ds} -subgroup moderna -out epitope_wt_thresholds-moderna.${params.plot_format}
        """ 
}

haa = Channel.fromPath("../analysis-scripts/haarvi-subgroups.py").combine(bde[3])

process haarvi_subgroup {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(haarvi),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from haa

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${haarvi} -dataset ${layered_phip_ds} -batch ${batch} -out haarvi.${params.plot_format}
        """ 
}

nihc = Channel.fromPath("../analysis-scripts/nih-subgroup.py").combine(bde[4])

process nih_subgroup {
    
    publishDir "${params.phip_data_dir}/${params.alignment_tool}/${batch}/${params.plot_format}", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 
    

    input:
        set(
            file(nih),
            val(batch),
            file(layered_phip_ds),
            file(epitopes)
        ) from nihc

    output:
        file "*.${params.plot_format}"

    script:
        """
        set -eu
        python ${nih} -dataset ${layered_phip_ds} -batch ${batch} -out nih.${params.plot_format}
        """ 
}

