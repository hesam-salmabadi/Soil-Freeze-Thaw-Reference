"""CF-compliant NetCDF writer for the release dataset.

Serializes per-site time series (soil temperature, effective permittivity,
freezing probability, soil-state classification) plus fitted/pooled SFCC
parameters and their uncertainties into a CF-compliant NetCDF file with
appropriate coordinates, units, and metadata attributes.

TODO: implement ``write_netcdf(dataset, path)`` (xarray -> netCDF4).
"""
