#!/usr/bin/env python
# coding: utf-8

import numpy as np, camb, basic, tqdm
from scipy.interpolate import interp1d

# Function to get the sigma_8 and P(k) at z = 0.0 from CAMB
def get_pklin0(h0,Ob,Om,As,ns,maxkh=10.,k_per_logint=200):

    # Reset CAMB
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=1e2*h0, ombh2=Ob*h0**2, omch2=(Om-Ob)*h0**2,omk=0,mnu=0)
    pars.set_dark_energy() #LCDM (default)
    pars.InitPower.set_params(ns=ns, r=0, As=As)

    # Set matter power spectrum at z=0: P(k,z)
    pars.set_matter_power(redshifts=[0.],kmax=maxkh/h0*1.1,k_per_logint=k_per_logint)
    pars.NonLinear = camb.model.NonLinear_none
    
    # Calculate the intermediate results for these parameters
    results = camb.get_results(pars)
    sigma8 = results.get_sigma8()[0]
    
    #print('sigma8',sigma8)
    
    # Calculate the CAMB power spectra to be interpolated
    kh, _z, pk = results.get_linear_matter_power_spectrum(have_power_spectra=True,nonlinear=False)

    return kh, pk[0], sigma8


def Lpoints(lmin,lmax,bn,xmin=-np.pi/2.,xmax=np.pi/2.):
    x = np.arange(xmin,xmax+(xmax-xmin)/bn,(xmax-xmin)/bn)
    b = np.sin(x)
    db = 0.*b
    db[1:] = (b[1:]-b[:-1])*(lmax-lmin)/(b[-1]-b[0])
    l = lmin + np.cumsum(db)
    l = np.unique(l.astype(int))
    if l[-1]!=lmax: l[-1] = lmax
    return l


def Lpoints2(lmin,lmax,bn,mid0=400,mid1=1600):
    bn0 = np.int(bn*1./2.)
    bn2 = np.int(bn*1./3.)
    bn1 = bn-bn0-bn2
    print(bn0,bn1,bn2)
    l0 = Lpoints(lmin,mid0,bn0,xmax=0.)
    l1 = Lpoints(mid0,mid1,bn1,xmax=0.)
    l2 = Lpoints(mid1,lmax,bn2,xmin=0.)
    l  = np.unique(np.concatenate((l0,l1,l2)))
    return l


def compute_skewspec(i,j,D='data_local/skewspec/',zn=20,zss=1.0334,lmin=100,lmax=2000,bn=40,Om=.3,h0=.7,Ob=0.05,As=2e-9,ns=0.96,verbose=False):

    # interpolation points
    if bn is None:
        bn = lmax-lmin
        oL = np.linspace(lmin,lmax,bn+1,dtype=np.int)
    else:
        #ol = Lpoints(lmin,lmax,bn)
        oL = Lpoints2(lmin,lmax,bn)
    # output multipoles
    L  = np.linspace(lmin,lmax,lmax-lmin+1)
    
    # Pk
    kh, pklin0, s8 = get_pklin0(h0,Ob,Om,As,ns)

    # Sl
    zmin, zmax = 0.0001, min(zss,zss)
    skew = basic.bispec.skewspeclens('input','RT',zmin,zmax,zn,[zss,zss],oL,lmin,lmax,kh,pklin0,pb=True,Om=Om,H0=h0*1e2,mnu=0.,ns=ns,verbose=verbose)

    # interpolate
    sl = np.zeros((3,lmax-lmin+1))
    for s in range(3):
        sl[s] = interp1d(oL, skew[s], kind='cubic')(L)
    
    # results
    cp = {'h0':h0,'Omega_m':Om,'Omega_b':Ob,'As':As,'ns':ns,'sigma8':s8}
    Sl = {'L':L,'Sl0':sl[0],'Sl1':sl[1],'Sl2':sl[2]}
    
    f = D + 'Sl_grid'+str(i)+'-'+str(j)+'_zs'+str(zss)+'_zn'+str(zn)+'_l'+str(lmin)+'-'+str(lmax)+'_'+str(bn)
    np.savez(f,cp=cp,Sl=Sl)


def params_array(pmin,pmax,pnum):
    dp = (pmax-pmin)/pnum
    return np.arange(pmin,pmax+dp,dp)
    

# S8 = s8 x (Om/0.3)**0.5 = As**(1/2) * Om**(5/4) * h**(7/4)
# Mh = (As/As0) * (Om/0.3)**2
Ommin, Ommax, Omnum = .2, .4, 50
Mhmin, Mhmax, Mhnum = .5, 1.2, 50

for i, Om in enumerate(params_array(Ommin,Ommax,Omnum)):
    
    for j, Mh in enumerate(tqdm.tqdm(params_array(Mhmin,Mhmax,Mhnum))):
        
        #compute_skewspec(i,j,D='data_local/skewspec/test/',zn=30,lmin=100,lmax=2000,bn=40,Om=Om,As=Mh*(0.3/Om)**2*2e-9,verbose=True)
        compute_skewspec(i,j,D='data_local/skewspec/emu/',zn=30,lmin=100,lmax=2000,bn=40,Om=Om,As=Mh*(0.3/Om)**2*2e-9,verbose=True)

