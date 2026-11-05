"""Compatibility imports for the plotting application's equilibrium reader."""

from equilibrium_io import ASYM_FAM_NAMES, FAM_NAMES, load_state as _load_state


def load_state(path):
    try:
        return _load_state(path)
    except ValueError as error:
        raise SystemExit(str(error)) from error
