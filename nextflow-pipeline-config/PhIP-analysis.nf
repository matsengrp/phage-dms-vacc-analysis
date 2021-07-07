/*
PhIP-Flow Pipeline Script

Author: Jared G. Galloway
*/


peptide_reference_ch = Channel.fromPath("${params.peptide_table}")

// CONVERT PEPTIDE METADATA TO FASTA
process generate_fasta_reference {

    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/phippery:vacc-ms-analysis' 
    //publishDir "${params.phip_data_dir}/"

    input: file "pep_ref" from peptide_reference_ch

    output: file "peptides.fasta" into pep_channel_fasta

    shell:    
        """
        phippery peptide-md-to-fasta -d ${pep_ref} -o peptides.fasta
        """
}


// GENERATE INDEX
process generate_index {
 
    //publishDir "${params.phip_data_dir}/reference"
    container = 'quay.io/jgallowa/bowtie2:latest' 
    label 'multithread'

    input:
        file "pep_fasta" from pep_channel_fasta

    output:
        set(
            val("peptide_ref"), 
            file("peptide_index") 
        ) into pep_channel_index

    shell:    

        //if ("$params.alignment_tool" == "bowtie")
        //    """
        //    mkdir peptide_index
        //    bowtie-build --threads 4 \
        //    ${pep_fasta} peptide_index/peptide
        //    """

        //else if ("$params.alignment_tool" == "bowtie2")
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
    container = 'quay.io/jgallowa/bowtie2:latest' 

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

        """
        set -eu
        ${params.fastq_stream_func} ${respective_replicate_path} | \
        bowtie2 ${params.align_args} -x ${index}/peptide - > ${sample_id}.sam
        """
}


// SPLIT CHANNEL FOR COUNTS AND STATS IN PARALLEL
aligned_reads_sam.into{aligned_reads_for_counts; aligned_reads_for_stats}


// COMPUTE ALIGNMENT STATS FOR ALL STATS
process sam_to_stats {

    label 'multithread'
    container = 'quay.io/biocontainers/samtools:1.3--h0592bc0_3'

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
    container = 'quay.io/biocontainers/samtools:1.3--h0592bc0_3'

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
    
    publishDir "${params.phip_data_dir}/", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/phippery:vacc-ms-analysis' 

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
    
    publishDir "${params.phip_data_dir}/", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/phippery:latest' 

    input:
        file phip_ds from phip_data_ch
        //file alignment-stats from Channel.fromPath("../analysis-image/alignment-stats.py")
        file analysis from Channel.fromPath("../analysis-image/layer-stats.py")
        

    output:
        file "layered-analysis.phip" into layered_phip_data_ch

    """
    set -eu
    python ${analysis} ${phip_ds} layered-analysis.phip
    """ 
    //python compute_all_analysis_all_samples.py ${phip_ds} alignment-stats.pdf
}

// RUN ALL ANALYSIS
process analysis_plotting {
    
    publishDir "${params.phip_data_dir}/", mode: 'copy'
    label 'single_thread_large_mem'
    container = 'quay.io/matsengrp/vacc-ms-analysis:latest' 

    input:
        file layered_phip_ds from layered_phip_data_ch

    //output:
    //    file "${params.dataset_prefix}.phip" into phip_data_ch

    script:
        """
        set -eu
        python pca-scatter-directions.py ${layered_phip_ds} pca.pdf
        python heatmap_boxplot.py ${layered_phip_ds} heatmap-boxplot.pdf
        python logopairs-boxplot.py ${layered_phip_ds} logopairs-boxplot.pdf
        python haarvi-subgroups.py ${layered_phip_ds} haarvi.pdf
        python nih-subgroup.py ${layered_phip_ds} nih.pdf
        """ 
}
