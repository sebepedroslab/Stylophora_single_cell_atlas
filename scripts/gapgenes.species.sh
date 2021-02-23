#!/bin/bash
#$ -N gapgenes_species
#$ -t 1-13
#$ -cwd
#$ -q long-sl7,mem_512_12h 
#$ -M anamaria.elek@crg.eu
#$ -m a
#$ -o logs/gapgenes_species_$TASK_ID.out
#$ -j y
#$ -pe smp 12

wdir="/users/asebe/aelek/proj/scRNAseq_spis/Stylophora_single_cell_atlas/cross_species/tree"
tree=${wdir}"/species_vargenes1.8_coocaveragep0.75hclusth0.75,0.95_politomies0.10.RDS"
gene_expression=${wdir}"/Spis_Nvec_Xesp_Hvul_broad_cell_type_gene_FC"
outdir=${wdir}"/data"
arr=(0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2)
bg=${arr[${SGE_TASK_ID}]}

echo "ABS 0 ABSBG 0"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0 0
wait

echo "ABS 0 ABSBG 0.05"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0 0.05
wait

echo "ABS 0 ABSBG 0.1"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0 0.1
wait


echo "ABS 0.05 ABSBG 0"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.05 0
wait

echo "ABS 0.05 ABSBG 0.05"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.05 0.05
wait

echo "ABS 0.05 ABSBG 0.1"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.05 0.1
wait


echo "ABS 0.1 ABSBG 0"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.1 0
wait

echo "ABS 0.1 ABSBG 0.05"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.1 0.05
wait

echo "ABS 0.1 ABSBG 0.1"
Rscript --vanilla gapgenes.R $tree $gene_expression $outdir $bg 0.1 0.1
wait

echo "Done!"
