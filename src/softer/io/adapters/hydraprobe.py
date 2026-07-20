"""HydraProbe (Stevens Water) adapter — Sect. 2.1.1, Appendix A.

Networks: Candle Lake (BT), Kenaston (KN), Trail Valley Creek (TV).
Reports soil temperature and moisture together. Emits rows in the common
schema, preserving the probe's native output alongside soil temperature; the
raw-output -> effective-permittivity conversion is applied in
:mod:`softer.preprocess.permittivity`.

TODO: implement ``HydraProbeAdapter(BaseAdapter)``.
"""
