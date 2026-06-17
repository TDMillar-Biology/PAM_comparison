project="."

mkdir -p $project/workflow/rules
mkdir -p $project/workflow/envs
mkdir -p $project/workflow/scripts
mkdir -p $project/config
mkdir -p $project/data
mkdir -p $project/results
mkdir -p $project/logs
mkdir -p $project/containers


touch $project/workflow/Snakefile
touch $project/workflow/rules/mapping_qv.smk
touch $project/workflow/envs/mapping.yaml
touch $project/workflow/scripts/compute_qv.py
touch $project/config/config.yaml
touch $project/config/master_control.tsv
touch $project/README.md
touch $project/.gitignore

echo "Project scaffold created at $project/"
