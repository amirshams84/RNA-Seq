# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: July-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-seq
# Aim: Conda script to build environment for RNA-seq pipeline
# ################################### RECIPE ##################################### #

# ++++++++++++++++++++++++++++++
# Step1: CREATE BASIC CONDA ENVIRONMENT
# ------------------------------

conda update -n base -c defaults conda

conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

conda create -n RNA-seq_conda -c r python==3.7.0 r-base==3.5.1 -y

# ++++++++++++++++++++++++++++++
# Step2: Activate the environment
# ------------------------------

conda activate RNA-seq_conda

# ++++++++++++++++++++++++++++++
# Step3: Install tidyverse
# ------------------------------

conda install -c conda-forge r-tidyverse -y

# ++++++++++++++++++++++++++++++
# Step4: Install snakemake
# ------------------------------

conda install -c bioconda snakemake -y

# ++++++++++++++++++++++++++++++
# Step5: Install applications
# ------------------------------

conda install -c bioconda samtools bamtools star=2.7.0f star-fusion=1.6.0 fastp fastqc multiqc bedtools deeptools mosdepth subread picard bcftools minimap2 kraken2 krona rseqc htseq bioconductor-dupradar -y

conda install -c conda-forge r-kernsmooth -y

conda install -c anaconda gfortran_impl_linux-64==7.5.0 gxx_impl_linux-64==7.5.0 gcc_impl_linux-64==7.5.0 -y

# ++++++++++++++++++++++++++++++
# Step6: Export
# ------------------------------
conda env export --no-builds | grep -v "^prefix: " > RNA-Seq_conda_nobuild.yaml

