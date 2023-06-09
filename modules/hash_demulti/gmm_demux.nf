#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process gmm_demux{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/gmm_demux", mode:'copy'
    input:
        each path_hto
        //HTO names as string separated by commas
        each hto_name_gmm
        //5 cases are available for the tool it all depends on the vars given
        //mode 2
        //need estimate number of cells in the single cell assay
        //obligatory
        each summary
        //need to be combined with summary to get a report as file
        each report 
        //mode 4
        // write csv or tsv - type of input
        each mode
        //case 5
        each extract 
        //float between 0 and 1
        each threshold_gmm
        
        
    
    output:
        path "gmm_demux_${task.index}"
        
    script:
        def extract_droplets = extract != 'None' ? " -x ${extract}" : ''
        
        if(mode=="csv")
            """
            mkdir gmm_demux_${task.index}
            
            GMM-demux -c $path_hto $hto_name_gmm -u $summary --report demuxem_${task.index} --full demuxem_${task.index} $extract_droplets -t $threshold_gmm
            """
        else
            """
            mkdir gmm_demux_${task.index}
            
            GMM-demux $path_hto $hto_name_gmm -u $summary -r gmm_demux_${task.index}/$report --full gmm_demux_${task.index} -o gmm_demux_${task.index} $extract_droplets -t $threshold_gmm
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

workflow gmm_demux_hashing{
  main:
        path_hto = split_input(params.hto_matrix_preprocess)
        hto_name_gmm = split_input(params.hto_name_gmm)
        summary = split_input(params.summary)
        report = split_input(params.report)
        mode = split_input(params.mode)
        extract = split_input(params.extract)
        threshold_gmm = split_input(params.threshold_gmm)

        gmm_demux(path_hto,hto_name_gmm,summary,report,mode,extract,threshold_gmm)
  
  emit:
        gmm_demux.out.collect()
}


workflow{
    gmm_demux_hashing()

}
