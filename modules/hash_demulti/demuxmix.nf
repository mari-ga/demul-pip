#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process demuxmix{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/demuxmix", mode:'copy'
    input:
        //shares pre-process Seurat
        each seurat_object
        //Same assay as Seurat
        each assay
        each model
        each alpha_demuxmix
        each beta_demuxmix
        each tol_demuxmix
        each maxIter_demuxmix
        each k_hto
        each k_rna
        
    output:
        path "demuxmix_${task.index}"
        
    script:
        def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot ${generate_gender_plot}" : ''
        """
        mkdir demuxmix_${task.index}
        demuxmix.R --seuratObject $seurat_object --assay $assay --model $model --alpha_demuxmix $alpha_demuxmix \
            --beta_demuxmix $beta_demuxmix --tol_demuxmix $tol_demuxmix --maxIter_demuxmix $maxIter_demuxmix \
            --k_hto $k_hto --k_rna $k_rna --outputdir demuxmix_${task.index}
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

workflow demuxmix_hashing{
  take:
        seurat_object
  main:
        assay = split_input(params.assay)
        model = split_input(params.model)
        alpha_demuxmix =  split_input(params.alpha_demuxmix)
        beta_demuxmix = split_input(params.beta_demuxmix)
        tol_demuxmix = split_input(params.tol_demuxmix)
        maxIter_demuxmix = split_input(params.maxIter_demuxmix)
        k_hto = split_input(params.k_hto)
        k_rna = split_input(params.k_rna) 

        demuxmix(seurat_object,assay,model, alpha_demuxmix, beta_demuxmix, tol_demuxmix, maxIter_demuxmix, k_hto, k_rna )
  
  emit:
        demuxmix.out.collect()
}

workflow{
    demuxmix_hashing()
}
