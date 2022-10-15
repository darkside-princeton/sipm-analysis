from scipy.special import erf, gamma
import numpy as np
import scipy.integrate as integrate
import sipm.constants as const

NORM_SR90 = 6004718046792.3545
NORM_Y90 = 2.056512394298567e+16

def Fermi_function(T: float, Z: int, A: int) -> float:
    '''s
    Z: atomic number of final state nucleus
    A: mass number of final state nucleus
    T: kinetic energy in keV
    '''
    E = const.M_ELE + T #keV
    p = np.sqrt(E**2 - const.M_ELE**2) #keV
    S = np.sqrt(1-(const.ALPHA*Z)**2)
    eta = const.ALPHA*Z*E/p
    R_N = A**(1/3)*1.2e-15 #m
    return 2*(1+S)/gamma(1+2*S)**2 * (2*p*R_N/const.HBAR/const.C)**(2*(S-1)) * np.exp(np.pi*eta) * abs(gamma(complex(S, eta)))**2

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
        E = const.M_ELE + T #keV
        p = np.sqrt(E**2-const.M_ELE**2) #keV
        if u1f:
            return C*Fermi_function(E,Z,A)*p*E*(Q-T)**2*U1F_shape(W=1+T/const.M_ELE, W0=1+Q/const.M_ELE)
        else:
            return C*Fermi_function(E,Z,A)*p*E*(Q-T)**2

def beta_shape_norm(Q: float, Z: int, A: int, u1f: bool) -> float:
    return integrate.quad(lambda x: beta_shape(x, Q, 1, Z, A, u1f), 0, Q)

def response(E: float, Npe: float, gamma: float, Lyp: float) -> float:
    return np.exp(-0.5*(Npe-E*Lyp)**2/(Npe*gamma))/np.sqrt(2*np.pi*Npe*gamma)

def beta_spectrum(Npe: float, C: float, gamma: float, Lyp: float, Q: float, Z: int, A: int, u1f: bool) -> float:
    return np.array([integrate.quad(lambda x: beta_shape(x, Q, C, Z, A, u1f)*response(x, Npe_, gamma, Lyp), 0, Q)[0] for Npe_ in Npe])

def trigger_eff(Npe: float, N0: float, B: float) -> float:
    '''
    N0: cutoff
    B: sharpness
    '''
    return 0.5*(1+erf((Npe-N0)/B/np.sqrt(Npe)))

def backgrounds(Npe: float, b=float) -> float:
    '''
    b=slope
    '''
    return b*np.exp(-b*Npe)

def model(Npe: float, n_beta:float, n_bkg:float, N0_trig:float, B_trig:float, enf:float, b_bkg:float, lyp: float):
    '''
    Count [1/Npe]
    n_beta: beta count
    n_bkg: background count
    N0_trig, B_trig: trigger efficiency parameters
    b_bkg: background parameters
    enf: excess noise factor
    lyp: gross light yield 
    return: sr90, y90, bkg, total
    '''
    sr90 = 0.5*beta_spectrum(Npe, C=1/NORM_SR90, gamma=enf, Lyp=lyp, Q=const.Q_SR90, Z=const.Z_Y90, A=const.A_SR90, u1f=True)
    y90 = 0.5*beta_spectrum(Npe, C=1/NORM_Y90, gamma=enf, Lyp=lyp, Q=const.Q_Y90, Z=const.Z_ZR90, A=const.A_SR90, u1f=True)
    bkg = backgrounds(Npe, b=b_bkg)
    trig = trigger_eff(Npe, N0=N0_trig, B=B_trig)
    return trig*n_beta*sr90, trig*n_beta*y90, trig*n_bkg*bkg, trig*(n_beta*(sr90+y90)+n_bkg*bkg)

def model_total(Npe: float, n_beta:float, n_bkg:float, N0_trig:float, B_trig:float, enf:float, b_bkg:float, lyp: float):
    '''
    Count [1/Npe]
    n_beta: beta count
    n_bkg: background count
    N0_trig, B_trig: trigger efficiency parameters
    b_bkg: background parameters
    enf: excess noise factor
    lyp: gross light yield 
    return: sr90, y90, bkg, total
    '''
    sr90 = 0.5*beta_spectrum(Npe, C=1/NORM_SR90, gamma=enf, Lyp=lyp, Q=const.Q_SR90, Z=const.Z_Y90, A=const.A_SR90, u1f=True)
    y90 = 0.5*beta_spectrum(Npe, C=1/NORM_Y90, gamma=enf, Lyp=lyp, Q=const.Q_Y90, Z=const.Z_ZR90, A=const.A_SR90, u1f=True)
    bkg = backgrounds(Npe, b=b_bkg)
    trig = trigger_eff(Npe, N0=N0_trig, B=B_trig)
    return trig*(n_beta*(sr90+y90)+n_bkg*bkg)
