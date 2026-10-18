# 激活 openmm_dev（MVP-0 计划 §1.1）
# 用法: source scripts/env.sh
export MAMBA_EXE="${HOME}/miniforge3/bin/mamba"
export MAMBA_ROOT_PREFIX="${HOME}/miniforge3"
source "${HOME}/miniforge3/etc/profile.d/mamba.sh"
mamba activate openmm_dev

# Load the user's existing engine module in this shell only.
if [ -r "${HOME}/modulefiles/gromacs2026.3" ]; then
    source /etc/profile.d/env-modules.sh
    module use "${HOME}/modulefiles"
    module load gromacs2026.3
fi
