# Trimal as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trimal-fasttree --aligner clustalo --trimmer trimal --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trimal-phyml --aligner clustalo --trimmer trimal --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trimal-raxml --aligner clustalo --trimmer trimal --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trimal-iqtree --aligner clustalo --trimmer trimal --tree_builder iqtree -work-dir result/work/

# Clipkit as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-clipkit-fasttree --aligner clustalo --trimmer clipkit --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-clipkit-phyml --aligner clustalo --trimmer clipkit --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-clipkit-raxml --aligner clustalo --trimmer clipkit --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-clipkit-iqtree --aligner clustalo --trimmer clipkit --tree_builder iqtree -work-dir result/work/

# trim_alg_v2 as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trim_alg_v2-fasttree --aligner clustalo --trimmer trim_alg_v2 --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trim_alg_v2-phyml --aligner clustalo --trimmer trim_alg_v2 --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trim_alg_v2-raxml --aligner clustalo --trimmer trim_alg_v2 --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/clustalo-trim_alg_v2-iqtree --aligner clustalo --trimmer trim_alg_v2 --tree_builder iqtree -work-dir result/work/
