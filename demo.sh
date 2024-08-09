nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/ --aligner mafft --trimmer trimal --tree_builder fasttree -work-dir result/work/ --customConfig workflow.config
