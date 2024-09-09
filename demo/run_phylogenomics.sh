#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8GB
#SBATCH --array=1-18495
#SBATCH -e /scratch/plaza/projects/EggNOG6_PfamFams_v35/slurm_err/phylogenomics/pfamfam_phylogenomics_%A_%a.err
#SBATCH -o /scratch/plaza/projects/EggNOG6_PfamFams_v35/slurm_out/phylogenomics/pfamfam_phylogenomics_%A_%a.out


#conda activate  eggnog6_phylogenetics

export PATH=/scratch/plaza/soft/FAMSA_v2:$PATH
export PATH=/scratch/plaza/soft:$PATH
raw_fasta=$(cat /scratch/plaza/projects/EggNOG6_PfamFams_v35/path2fastas.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

name=$(echo $raw_fasta | cut -d '/' -f 7)

num_seqs=$(grep -c '^>' $raw_fasta)

echo $name $num_seqs

start=`date +%s`
if [ "$num_seqs" -gt 1000 ]; then
	echo famsa
	famsa $raw_fasta /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/alg/$name.aln

else
	echo mafft
	mafft --auto --anysymbol  $raw_fasta > /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/alg/$name.aln

fi
end_1=`date +%s`


python /scratch/plaza/projects/EggNOG6_PfamFams_v35/scripts/trim_alg_v2.py -i /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/alg/$name.aln -o /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/trim_alg/$name.trim.aln --min_res_percent 0.1  --min_res_abs 3
end_2=`date +%s`

FastTree /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/trim_alg/$name.trim.aln > /scratch/plaza/projects/EggNOG6_PfamFams_v35/phylogenomics/trees/$name.nw
end_3=`date +%s`
echo Execution_1:`expr $end_1 - $start`
echo Execution_2:`expr $end_2 - $end_1`
echo Execution_3:`expr $end_3 - $end_2`
echo Total_Execution:`expr $end_3 - $start`