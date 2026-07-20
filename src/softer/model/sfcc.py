"""SFCC model in permittivity-temperature space (Sect. 2.3.1, Eq. 4).

Recasts the Amankwah et al. liquid-water-content SFCC (Eq. 1) into effective
permittivity via a dielectric mixing model (Eqs. 2-3), giving eps_eff as a
function of soil temperature T:

    eps_eff(T) = ( (eps_int**a - eps_res**a) * exp(b * (T - T_f)) + eps_res**a )**(1/a)   if T <  T_f
    eps_eff(T) =   eps_int                                                                  if T >= T_f

Parameters:
    - eps_int : pre-freezing effective permittivity (total/initial water content)
    - eps_res : residual effective permittivity (residual water content)
    - b       : transition sharpness / shape factor [degC^-1]
    - T_f     : freezing-onset temperature [degC]
    - a (alpha): dielectric mixing exponent, fixed at 0.5 (low sensitivity;
                 see Sect. 2.3.2).

Because state classification (frozen/transitional/unfrozen) is qualitative, the
permittivity -> water-content conversion is unnecessary here.

TODO: implement ``sfcc(T, eps_int, eps_res, b, T_f, alpha=0.5)`` (vectorized,
piecewise) plus a ``normalized`` variant used for freezing probability.
"""
