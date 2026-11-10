# 激活 openmm_dev（MVP-0 计划 §1.1）
# 用法: source scripts/env.sh
# 可覆盖（决策 0026）：FFC_MAMBA_ROOT（默认 ~/miniforge3）、FFC_ENV_NAME（默认 openmm_dev）、
#                     FFC_MODULEFILES（默认 ~/modulefiles）、FFC_GMX_MODULE（默认 gromacs2026.3）
_ffc_root="${FFC_MAMBA_ROOT:-${HOME}/miniforge3}"
export MAMBA_EXE="${_ffc_root}/bin/mamba"
export MAMBA_ROOT_PREFIX="${_ffc_root}"
source "${_ffc_root}/etc/profile.d/mamba.sh"
mamba activate "${FFC_ENV_NAME:-openmm_dev}"

# Load the user's existing engine module in this shell only.
_ffc_mod_dir="${FFC_MODULEFILES:-${HOME}/modulefiles}"
_ffc_mod="${FFC_GMX_MODULE:-gromacs2026.3}"
if [ -r "${_ffc_mod_dir}/${_ffc_mod}" ]; then
    source /etc/profile.d/env-modules.sh
    module use "${_ffc_mod_dir}"
    module load "${_ffc_mod}"
fi
unset _ffc_root _ffc_mod_dir _ffc_mod
