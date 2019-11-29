# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""
from astropy.io import fits
import numpy as np 
import math
import sys, os
import matplotlib.pyplot as plt



def geocutoff_cumani(h):                         #cutoff a latitudine costante
    l_const = 0*np.pi/180                        #latitudine 0° in radianti
    R_earth = 6371                               #km
    return ((11.9*((R_earth+550)/(R_earth+h))**2)*np.cos(l_const)**4)

def geocutoff_cumani2(l):
    l_rad = l*np.pi/180
    return (11.9*np.cos(l_rad)**4)               #equazione assottigliata ((num/den) = 1)
    


    
    
h = np.linspace(400,600,20)               #ascisse altitudine
l = np.linspace(0,50,20)                  #ascisse latitudine
print len(h)
print len(l)

plt.plot(h,geocutoff_cumani(h))
plt.show()
plt.plot(l,geocutoff_cumani2(l))
plt.show()




"""    
#def geocutoff_l(l):                      #funzione x geocutoff a h costante
#    h_const = 550*10**3                  #550km (in metri)
#    g10 = 29442.0*10**(-9)               #T, aggiornato al 2015
#    R = 6371*10**3                       #raggio Terra (m)
#    g1 = (g10*R)/4
#    g2 = (1+(h_const/R))**(-2)
#    g3 = (np.cos(l*np.pi/180)**4)
#    geocut_l = g1*g2*g3
#    return geocut_l

#def geocutoff_h(h):                     #funzione x geocutoff a l costante
#    l_const = 0*np.pi/180               #0°
#    g10 = 29442.0*10**(-9)              
#    R = 6371*10**3
#    g1 = (g10*R)/4
#    g2 = (1+(h/R))**(-2)
#    g3 = (np.cos(l_const)**4)
#    geocut_h = g1*g2*g3                 #inserito lo stesso g3 anche se = 1 per l=0° x eventuali modifiche future
#    return geocut_h

#def unmod_proton_flux(Ek):
#    Ek_converter = Ek*10**9            #converte i GeV in eV
#    E0 = 938*10**6                     #energia a riposo del protone
#    A = 23.9*10**(-6)                  #conteggi (s^-1 m^-2 sr^-1 eV^-1)
#    a = 2.83                           #adimensionale
#    AZ = 1                             #nucleoni/numero atomico
#    e = 1.602*10**(-19)                #carica fondamentale elettrone
#    Rigidity = ((AZ/e)**2)*(Ek_converter*(Ek_converter + 2 * E0))**(0.5)
#    unmod_Pflux = A * (Rigidity)**(-a)
#    return unmod_Pflux

#def solarmod_proton_flux(Ek):
#    Ek_converter = Ek*10**9            #converte i GeV in eV
#    phi = 800*10**6                    #modulazione solare (V)
#    M = 9.3828*10**8                   #massa protone (eV/c^2)
#    E0 = 938*10**6                     #energia a riposo del protone
#    A = 23.9*10**(-6)                  #conteggi (s^-1 m^-2 sr^-1 eV^-1)
#    a = 2.83                           #adimensionale
#    AZ = 1                             #nucleoni/numero atomico
#    e = 1.602*10**(-19)                #carica fondamentale elettrone
#    c = 299792458
#    Mc = M*(c)**2
#    Rigidity = ((AZ/e)**2)*(Ek_converter*(Ek_converter + 2 * E0))**(0.5)
#    unmod_Pflux = A * (Rigidity)**(-a)
#    solar_modulation = ((Ek_converter + Mc)**2 - (Mc)**2)/(((Ek_converter + Mc + e*phi)**2)-(Mc)**2)
#    solarmod_Pflux = unmod_Pflux*solar_modulation
#    return solarmod_Pflux

#def modulated_proton_flux(Ek):
#    Ek_converter = Ek*10**9            #converte i GeV in eV
#    phi = 800*10**6                    #modulazione solare (V)
#    h = 38000                          #altitudine (m)
#    l = 0.735                          #latitudine geomagnetica (rad)
#    M = 9.3828*10**8                   #massa protone (eV/c^2)
#    E0 = 938*10**6                     #energia a riposo del protone
#    A = 23.9*10**(-6)                  #conteggi (s^-1 m^-2 sr^-1 eV^-1)
#    a = 2.83                           #adimensionale
#    r = 12.0
#    AZ = 1                             #nucleoni/numero atomico
#    e = 1.602*10**(-19)                #carica fondamentale elettrone
#    c = 299792458
#    Mc = M*(c)**2
#    Rigidity = ((AZ/e)**2)*(Ek_converter*(Ek_converter + 2 * E0))**(0.5)
#    R_cutoff = 4.46*10**9
#    unmod_Pflux = A * (Rigidity)**(-a)
#    solar_modulation = ((Ek_converter + Mc)**2 - (Mc)**2)/(((Ek_converter + Mc + e*phi)**2)-(Mc)**2)
#    modulated_Pflux = solar_modulation * (1/(1+(Rigidity/R_cutoff))**(-r))
#    return modulated_Pflux

for i in [h]:                          #grafico zoom geocutoff a l costante
    plt.plot(h, geocutoff_h(h))
    plt.show()
    
#for i in [Ek]:
#    plt.plot(Ek, unmod_proton_flux(Ek))
#    plt.show
    
#for i in [Ek]:
#    plt.plot(Ek, solarmod_proton_flux(Ek))
#    plt.show
    
    
#Ek_converter = 1000
#phi = 800*10**6                    #modulazione solare (V)
#M = 9.3828*10**8                   #massa protone (eV/c^2)
#E0 = 938*10**6                     #energia a riposo del protone
#A = 23.9*10**(-6)                  #conteggi (s^-1 m^-2 sr^-1 eV^-1)
#a = 2.83                           #adimensionale
#AZ = 1                             #nucleoni/numero atomico
#e = 1.602*10**(-19)                #carica fondamentale elettrone
#c = 299792458
#Mc = M*(c)**2
#Rigidity = ((AZ/e)**2)*(Ek_converter*(Ek_converter + 2 * E0))**(0.5)
#unmod_Pflux = A * (Rigidity)**(-a)
#solar_modulation = (Ek_converter + Mc)**2 - (Mc)**2
#aaa = (Mc**2)**2
#ciao = 1
#print (Mc)
#print (Ek_converter)
#print (solar_modulation)
#print (aaa)
#print (ciao)

#Ek = np.linspace(0.001,100)            #range energie cinetiche

#for i in [l]:                          #grafico generale
#    plt.plot(l, geocutoff_l(l), geocutoff_h(h))
#    plt.show()

"""
    