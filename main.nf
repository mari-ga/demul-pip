#!/usr/bin/env nextflow
nextflow.enable.dsl=2
'''
BAM Index should in the format BAMFILE.bam.bai

'''
include { gene_demultiplexing } from './modules/gene_demultiplexing'
include { hash_demultiplexing } from './modules/hash_demultiplexing'

process summary_all{
    publishDir "$params.outdir/$params.mode/summary", mode: 'copy'
    input:
        path gene_demulti_result
        path hash_demulti_result
    output:
        path '*.csv'

    script:
        """
        summary.R --gene_demulti $gene_demulti_result --hash_demulti $hash_demulti_result
        """
}

process filter_barcodes{
    publishDir "$params.outdir/$params.mode/compare", mode: 'copy'
    input:
        each selected_barcodes
        each white_list
    output:
        path '*.tsv'

    script:
        
        """
        barcode_name="\$(basename \$(dirname $white_list))/\$(basename $white_list)"
        barcode_name="\${barcode_name////.}"
        comm -12 <(sort $selected_barcodes) <(sort $white_list) > \${barcode_name}
        """
}


def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}

workflow{
    mode = params.mode
    whitelist_barcodes = split_input(params.barcodes)
    if (mode == "genetic"){
        gene_demultiplexing(whitelist_barcodes)
    }
    else if (mode == "hash"){
        hash_demultiplexing()
    }
    else if (mode == "parallel"){
        hash_demultiplexing()
        gene_demultiplexing(whitelist_barcodes)
        gene_summary = gene_demultiplexing.out
        hash_summary = hash_demultiplexing.out
        summary_all(gene_summary, hash_summary)
        
    }
    else{
        hash_demultiplexing()
        selected_barcodes = hash_demultiplexing.out.map{ return it + "/selected_barcodes.tsv"}
        filter_barcodes(selected_barcodes, whitelist_barcodes)
        recover_barcodes= filter_barcodes.out.collect()
        gene_demultiplexing(recover_barcodes)
        gene_summary = gene_demultiplexing.out
        hash_summary = hash_demultiplexing.out
        summary_all(gene_summary, hash_summary)
    }
}
