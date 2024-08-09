# Trimal as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trimal-fasttree --aligner mafft --trimmer trimal --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trimal-phyml --aligner mafft --trimmer trimal --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trimal-raxml --aligner mafft --trimmer trimal --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trimal-iqtree --aligner mafft --trimmer trimal --tree_builder iqtree -work-dir result/work/

# Clipkit as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-clipkit-fasttree --aligner mafft --trimmer clipkit --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-clipkit-phyml --aligner mafft --trimmer clipkit --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-clipkit-raxml --aligner mafft --trimmer clipkit --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-clipkit-iqtree --aligner mafft --trimmer clipkit --tree_builder iqtree -work-dir result/work/

# trim_alg_v2 as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trim_alg_v2-fasttree --aligner mafft --trimmer trim_alg_v2 --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trim_alg_v2-phyml --aligner mafft --trimmer trim_alg_v2 --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trim_alg_v2-raxml --aligner mafft --trimmer trim_alg_v2 --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/mafft-trim_alg_v2-iqtree --aligner mafft --trimmer trim_alg_v2 --tree_builder iqtree -work-dir result/work/
