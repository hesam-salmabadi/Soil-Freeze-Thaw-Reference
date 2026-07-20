"""Command-line interface for SoFTeR.

Exposes the ``softer`` console script (see ``[project.scripts]`` in
pyproject.toml). The primary command is::

    softer run --config configs/config.example.yaml

which loads a run configuration and drives the end-to-end pipeline
(:mod:`softer.pipeline`).

TODO: implement the click command group. Suggested subcommands:
    - ``run``      end-to-end pipeline from a config file
    - ``ingest``   run only the ingest/adapter stage
    - ``fit``      run only the SFCC fitting stage
    - ``package``  build the CF-NetCDF release from processed outputs
"""
