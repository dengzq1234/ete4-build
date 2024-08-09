# Trimal as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trimal-fasttree --aligner famsa --trimmer trimal --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trimal-phyml --aligner famsa --trimmer trimal --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trimal-raxml --aligner famsa --trimmer trimal --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trimal-iqtree --aligner famsa --trimmer trimal --tree_builder iqtree -work-dir result/work/

# Clipkit as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-clipkit-fasttree --aligner famsa --trimmer clipkit --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-clipkit-phyml --aligner famsa --trimmer clipkit --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-clipkit-raxml --aligner famsa --trimmer clipkit --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-clipkit-iqtree --aligner famsa --trimmer clipkit --tree_builder iqtree -work-dir result/work/

# trim_alg_v2 as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trim_alg_v2-fasttree --aligner famsa --trimmer trim_alg_v2 --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trim_alg_v2-phyml --aligner famsa --trimmer trim_alg_v2 --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trim_alg_v2-raxml --aligner famsa --trimmer trim_alg_v2 --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/famsa-trim_alg_v2-iqtree --aligner famsa --trimmer trim_alg_v2 --tree_builder iqtree -work-dir result/work/
