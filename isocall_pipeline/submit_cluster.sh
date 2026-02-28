#!/bin/bash
#SBATCH --job-name=isocall_nf
#SBATCH --output=isocall_nf_%j.log
#SBATCH --error=isocall_nf_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -p linux12hr

# Nextflow head process — submits individual jobs via SLURM
# Run from: nmd/isocall_pipeline/

module load nextflow/24.10.4
module load java/21

PROJ_ROOT="/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025"

nextflow run main.nf \
    -profile slurm \
    --samplesheet assets/samplesheet.csv \
    --genome_fasta "${PROJ_ROOT}/results/reference_files/GRCh38.primary_assembly.genome.fa" \
    --annotation_gtf "${PROJ_ROOT}/results/reference_files/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf" \
    --outdir results \
    --output_prefix nmd_isocall \
    -resume
