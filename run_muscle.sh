# Trimal as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trimal-fasttree --aligner muscle --trimmer trimal --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trimal-phyml --aligner muscle --trimmer trimal --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trimal-raxml --aligner muscle --trimmer trimal --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trimal-iqtree --aligner muscle --trimmer trimal --tree_builder iqtree -work-dir result/work/

# Clipkit as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-clipkit-fasttree --aligner muscle --trimmer clipkit --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-clipkit-phyml --aligner muscle --trimmer clipkit --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-clipkit-raxml --aligner muscle --trimmer clipkit --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-clipkit-iqtree --aligner muscle --trimmer clipkit --tree_builder iqtree -work-dir result/work/

# trim_alg_v2 as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trim_alg_v2-fasttree --aligner muscle --trimmer trim_alg_v2 --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trim_alg_v2-phyml --aligner muscle --trimmer trim_alg_v2 --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trim_alg_v2-raxml --aligner muscle --trimmer trim_alg_v2 --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/muscle-trim_alg_v2-iqtree --aligner muscle --trimmer trim_alg_v2 --tree_builder iqtree -work-dir result/work/
