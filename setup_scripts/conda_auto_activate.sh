# System-wide Conda initialization and environment auto-activation
conda_base=`conda info --base`
CONDA_PROFILE="$conda_base/etc/profile.d/conda.sh"

if [ -f "$CONDA_PROFILE" ]; then
    . "$CONDA_PROFILE"
    # Only activate if it isn't already the active environment
    if [ "$CONDA_DEFAULT_ENV" != "snakemake" ]; then
        conda activate snakemake
    fi
fi
