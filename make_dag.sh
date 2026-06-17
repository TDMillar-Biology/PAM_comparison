snakemake --dag | dot -Tpng > dag.png
cat dag.tmp | dot -Tpng > figures/dag.png
