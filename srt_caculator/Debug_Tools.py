"""Debug utilities for inspecting solver parameters and runtime state."""

from __future__ import annotations

from typing import Mapping

from config import DEFAULT_PARAMS


def print_params(params: Mapping, title: str = "PARAMS") -> None:
    """Print a parameter mapping in a readable key-value format."""
    print(f"{title}:")
    for key, value in params.items():
        print(f"  {key}: {value}")


def print_default_params() -> None:
    """Print the default parameter list used by the solver."""
    print_params(DEFAULT_PARAMS, title="DEFAULT_PARAMS")
