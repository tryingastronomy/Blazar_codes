#!/usr/bin/env python2.7.16
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import csv

from matplotlib import rc, rcParams
rc('text',usetex=True)
#rc('axes', linewidth=2)
rc('font', weight='bold')
rc('font',**{'family':'serif','serif':['Computer Modern']})
rcParams['text.latex.preamble'] = r"\usepackage{boldmath} \usepackage{amsmath}"



def flux_from_magR(x):
    return 1e-23*3080.*10**(-x/2.5) #erg cm-2 s-1 hz-1

def lambda_to_nu(x):
    return 3e10/x

def normal(x, nu, index):
    return x/(nu**(index))

def alpha_ox(x, alpha):
    return x*10**(-alpha/0.384)



#For Z filter lambda = 913.4E-9m -> 913.4E-7cm
lambda_Z = 9.13E-5 #cm 

alpha = 1.4

lambda_R = 6.4e-5 #cm
lambda_ref = 2.4e-5 #cm
#nu_2kev = 2.42e17*2 #hz # if h = 4.135E-18 keV/Hz  so -> 1/h = 2.42E17 
nu_02kev = 2.42e17*3 #hz #change to XMM/NuSTAR etc. nu = E/h
nu_10kev = 2.42e17*79 #hz #change to XMM
#nu_03kev = 7.255E17




indexr=-0.5
indexx=-1.5


nu_R   = lambda_to_nu(lambda_Z)
nu_ref = lambda_to_nu(lambda_ref)

#name = np.loadtxt('list_of_source.txt', dtype=str, unpack=True, usecols=(0))
#magR = np.loadtxt('list_of_source.txt', unpack=True, usecols=(8))

#####
#-- For the ones with R mags
#####

magR = np.array([20.416])

fR = flux_from_magR(magR)

normR = normal(fR, nu_R, indexr)

f2500 = normR*(nu_ref**indexr)

f2kev = alpha_ox(f2500, alpha)

normx = normal(f2kev, nu_02kev, indexx)

f02to10kev = normx*integrate.quad(lambda x:x**(indexx), nu_02kev, nu_10kev)[0]

