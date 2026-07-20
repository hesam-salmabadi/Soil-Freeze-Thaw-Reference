"""TEROS12 (METER Group) adapter — Sect. 2.1.1, Appendix A.

Networks: James Bay (BJ), Montmorency Forest (FM), La Romaine (LR),
George River (GR). Reports soil temperature and moisture together. Emits rows in
the common schema; the raw-output -> effective-permittivity conversion is applied
in :mod:`softer.preprocess.permittivity`.

TODO: implement ``Teros12Adapter(BaseAdapter)``.
"""
