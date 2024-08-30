import numpy as np

def fourier_coefficient(eta, n, h_reflection):
    """
    Fourirer coefficients for rectangular signal.
    Needed for Barrier Bucket realization.

    Parameters
    ----------
    eta : float [1]
        Slippage factor. It is different for different particles
        and depent of 1) Energy [gamma]; 2) shift from reference [dp/p]
        eta = eta_0 + eta_1 * dp/p + ...
        The most interesting is 0, 1, 2 orders (tbh 0, 1)
    n : int [1]
        number of Fourirer coefficients
    h_reflection : float [1]
        Harmonic number of Barrier Bucket reflector.
        Usually defines as harmonic number of BB.
        defines as pi/phi_reflection,
        where phi_reflection = 2 * T_reflector / T_drift,
        T_reflector is duration of reflector signal [ns],
        T_drift is duration of drift space in BB [ns].

    For BB in NICA:
        eta change at transition at gamma_tr = 7.087,
        n = 1..25,
        T_reflector = 50 ns,
        T_drift = 400.25 ns. (Depend on f_revolution)
    """

    b = np.sign(eta)*2/(np.pi*n)*(1-np.cos((np.pi*n)/h_reflection))
    return b


def sigma_modulation(n, m, N):
    """
    Sigma_modulation procedure.
    [Unnormalized by PI! Used sinc(x / np.pi)
    to obtain the unnormalized sinc function]
    To make signal more common to the real.
    Using [sin(x)/x]^m.

    Parameters
    ----------
    n : int [1]
        number of Fourirer coefficient
    m : int [1]
        power of Sigma modulation
    N : int [1]
        total amount of Fourirer coefficient
    For BB in NICA:
        n = 1..25,
        m = 2 or 3
        N = 25
    """
    sigma = np.power(np.sinc((np.pi * n)/(2*(N+1))/np.pi), m)
    return sigma

def Voltage_amplitude(fourier_coefficient, sigma_modulation):
    Voltage_amplitude = sigma_modulation * fourier_coefficient
    return Voltage_amplitude
