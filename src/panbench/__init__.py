"""Stratified evaluation of variant calls against GIAB truth data.

The linear versus pangenome-aware caller comparison sketched in ``main.nf`` has not been
run. This package owns the evaluation, kept separate so the numbers can be tested without
a container runtime or a reference genome.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
