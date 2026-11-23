#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# config/init.sh
#
# Single entry point sourced by every script in this repository:
#
#   source "$(dirname "$0")/../config/init.sh"
#
# It loads the local configuration and the sample-selection helpers, and gives
# a readable error instead of a bare "No such file" when the repository has
# just been cloned.
# ---------------------------------------------------------------------------

_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "${_config_dir}/config.sh" ]]; then
    echo "ERROR: ${_config_dir}/config.sh not found." >&2
    echo "       Create it first:" >&2
    echo "         cp config/config.example.sh config/config.sh" >&2
    echo "       then edit the paths inside it." >&2
    exit 1
fi

source "${_config_dir}/config.sh"
source "${_config_dir}/select.sh"
