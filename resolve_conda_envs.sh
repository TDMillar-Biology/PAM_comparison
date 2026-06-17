# Compute nodes on HPRC do not have internet connection to resolve conda envs
# So we resolve them on the log in node before running the pipeline


module load GCC/12.3.0 OpenMPI/4.1.5 snakemake/8.4.2 Anaconda3/2024.02-1
snakemake --use-conda --conda-create-envs-only --conda-frontend conda


#source /sw/eb/sw/Anaconda3/2024.02-1/etc/profile.d/conda.sh

#conda create --prefix ./workflow/envs/svmu2 python=3.10

#conda run -p ./workflow/envs/svmu2 \
#    python -m pip install -e ./workflow/tools/svmu2

