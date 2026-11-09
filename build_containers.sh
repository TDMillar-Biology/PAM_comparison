for def in workflow/containers/definitions/*.def; do
    name=$(basename "$def" .def)

    sudo apptainer build \
        "workflow/containers/images/${name}.sif" \
        "$def" || exit 1
done
