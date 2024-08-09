# Trimal as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trimal-fasttree --aligner tcoffee --trimmer trimal --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trimal-phyml --aligner tcoffee --trimmer trimal --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trimal-raxml --aligner tcoffee --trimmer trimal --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trimal-iqtree --aligner tcoffee --trimmer trimal --tree_builder iqtree -work-dir result/work/

# Clipkit as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-clipkit-fasttree --aligner tcoffee --trimmer clipkit --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-clipkit-phyml --aligner tcoffee --trimmer clipkit --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-clipkit-raxml --aligner tcoffee --trimmer clipkit --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-clipkit-iqtree --aligner tcoffee --trimmer clipkit --tree_builder iqtree -work-dir result/work/

# trim_alg_v2 as trimmer
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trim_alg_v2-fasttree --aligner tcoffee --trimmer trim_alg_v2 --tree_builder fasttree -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trim_alg_v2-phyml --aligner tcoffee --trimmer trim_alg_v2 --tree_builder phyml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trim_alg_v2-raxml --aligner tcoffee --trimmer trim_alg_v2 --tree_builder raxml -work-dir result/work/
nextflow -log result/.nextflow.log run ete_build_dsl2.nf --input data/ --output result/tcoffee-trim_alg_v2-iqtree --aligner tcoffee --trimmer trim_alg_v2 --tree_builder iqtree -work-dir result/work/
