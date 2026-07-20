"""SoFTeR — Soil Freeze–Thaw Reference.

Automates the Soil Freezing Characteristic Curve (SFCC) framework of
Salmabadi et al. (2026, *The Cryosphere*, doi:10.5194/tc-20-1635-2026) and
applies it to heterogeneous in situ soil-temperature and soil-moisture-probe
records to produce a harmonized, non-binary soil freeze–thaw reference dataset.

The SFCC is fit in sensor-output-versus-temperature space (whatever quantity a
probe reports — permittivity, frequency, count, travel/output time, etc.),
not restricted to permittivity.
"""

__version__ = "0.0.1"
