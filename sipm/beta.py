from scipy.special import erf, gamma
import numpy as np
import scipy.integrate as integrate
Q_SR90 = 0.546e3 #keV
Q_Y90 = 2.28e3 #keV 
Z_Y90 = 39
Z_ZR90 = 40
A_SR90 = 90
ALPHA = 1/137
HBAR = 6.582e-19 #keV-s
C = 3e8 #m/s
M_ELE = 511 #keV

def Fermi_function(T: float, Z: int, A: int) -> float:
    '''s
    Z: atomic number of final state nucleus
    A: mass number of final state nucleus
    T: kinetic energy in keV
    '''
    E = M_ELE + T #keV
    p = np.sqrt(E**2 - M_ELE**2) #keV
    S = np.sqrt(1-(ALPHA*Z)**2)
    eta = ALPHA*Z*E/p
    R_N = A**(1/3)*1.2e-15 #m
    return 2*(1+S)/gamma(1+2*S)**2 * (2*p*R_N/HBAR/C)**(2*(S-1)) * np.exp(np.pi*eta) * abs(gamma(complex(S, eta)))**2

def U1F_shape(W: float, W0: float) -> float:
    '''
    ref: https://link.aps.org/doi/10.1103/PhysRev.82.48 Eq.1
    W: total energy in units of m_electron
    W0: maximum energy in units of m_electron
    '''
    return ((W0-W)**2+(W**2-1))/12

def U1F_shape_exact(W: float, W0: float, Z: int) -> float:
    '''
    ref: https://link.aps.org/doi/10.1103/PhysRev.82.48 Eq.2
    W: total energy in units of m_electron
    W0: maximum energy in units of m_electron
    Z: charge of final state nucleus
    TO DO
    '''
    return 0

def beta_shape(T: float, Q: float, C: float, Z: int, A: int, u1f: bool) -> float:
    '''
    T: kinetic energy in keV
    Q: maximum energy in keV
    C: normalization constant
    Z: atomic number of final state nucleus
    A: mass number of final state nucleus
    u1f: whether it's unique 1st forbidden decay or not
    '''
    if T>Q:
        return 0
    else:
        E = M_ELE + T #keV
        p = np.sqrt(E**2-M_ELE**2) #keV
        if u1f:
            return C*Fermi_function(E,Z,A)*p*E*(Q-T)**2*U1F_shape(W=1+T/M_ELE, W0=1+Q/M_ELE)
        else:
            return C*Fermi_function(E,Z,A)*p*E*(Q-T)**2

def beta_shape_norm(Q: float, Z: int, A: int, u1f: bool) -> float:
    return integrate.quad(lambda x: beta_shape(x, Q, 1, Z, A, u1f), 0, Q)

def response(E: float, Npe: float, gamma: float, Lyp: float) -> float:
    return np.exp(-0.5*(Npe-E*Lyp)**2/(Npe*gamma))/np.sqrt(2*np.pi*Npe*gamma)

def beta_spectrum(Npe: float, C: float, gamma: float, Lyp: float, Q: float, Z: int, A: int, u1f: bool) -> float:
    return integrate.quad(lambda x: beta_shape(x, Q, C, Z, A, u1f)*response(x, Npe, gamma, Lyp), 0, Q)[0]

def trigger_eff(Npe: float, N0: float, B: float) -> float:
    '''
    N0: cutoff
    B: sharpness
    '''
    return 0.5*(1+erf((Npe-N0)/B))
