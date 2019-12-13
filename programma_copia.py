# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""
from astropy.io import fits
from matplotlib.lines import Line2D
import numpy as np 
import math
import sys, os
import matplotlib.pyplot as plt


def geocutoff_cumani(h):                         # cutoff a latitudine costante
    l_const = 0*np.pi/180                        # latitudine 0° in radianti
    R_earth = 6371                               # km
    return ((11.9*((R_earth+550)/(R_earth+h))**2)*np.cos(l_const)**4)

def geocutoff_cumani2(l):
    l_rad = l*np.pi/180
    return (11.9*np.cos(l_rad)**4)               # equazione alleggerita ((num/den) = 1)

def g10(g):
    g10_values=[31543, 31464, 31354, 31212, 31060, 30926, 30805, 30715, 30654, 30594, 30554, 30500, 30421, 30334, 30220, 30100, 29992, 29873, 29775, 29692, 29619.4, 29554.63, 29496.57, 29442.0]
    return (g10_values)

def g10_var(g):
    g10_var=[79, 110, 142, 152, 134, 121, 90, 61, 60, 40, 54, 54, 79, 87, 114, 120, 108, 119, 98, 83, 72.6, 64.77, 58.06, 54.57]
    return (g10_var)
    
#def proton_flux(R):
#    C = 0.4544                                   # fattore di normalizzazione
#    gamma = -2.849                               # costante
#    dlt_gamma = 0.133                            # costante
#    s = 0.024                                    # smoothness of transition (?)
#    R0 = 336                                     # GV
#    A = 1                                        # numero di massa protone
#    Z = 1                                        # numero atomico protone
    #T0 = 938.25*10**(-3)                         # energia a riposo protone (GeV)
    #M0 = 1.6726219*10**(-27)                     # massa protone
    #arg = (T**2 + 2*T*T0)**(0.5)
    #arg_out = (T**2 + 2*T*M0)**(0.5)
    
    #R = (A/Z)*arg
    #R_out = arg_out = (T**2 + 2*T*M0)**(0.5)
    #F1 = (C*(R/45)**gamma)
    #F2 = (1+(R/R0)**(dlt_gamma/s))**s
    #F = ((C*(R/45)**gamma)*(1+(R/R0)**(dlt_gamma/s))**s)*R**2.7
    #f = (F1*F2)*R**(-2.7)
    
#    return (f)


   
h = np.linspace(400,600,20)                       # ascisse altitudine
l = np.linspace(0,50,20)                          # ascisse latitudine
g = np.linspace(1900,2015,24)                     # ascisse anno
#R = np.logspace(2,3,base=10.0)                    # ascisse energia cinetica (GeV)
#R = np.linspace(1,100000,10)

fig, ax1 = plt.subplots()
ax1.set_xlabel('Altitude / km', color = 'firebrick')
ax1.set_ylabel('Rcutoff / GV')
line1, = ax1.plot(h, geocutoff_cumani(h), '--', color = 'firebrick', label = 'Geomagnetic Latitude = 0$^{\circ}$')
ax1.tick_params(axis='x', labelcolor = 'firebrick')

ax2 = ax1.twiny()
color_blue = 'tab:blue'
ax2.set_xlabel('Geomagnetic Latitude / $^{\circ}$', color = color_blue)
line2, = ax2.plot(l, geocutoff_cumani2(l), color = color_blue, label = 'Altitude = 550 km')
ax2.tick_params(axis='x', labelcolor = color_blue)
ax1.grid(axis = 'y')
ax1.legend((line1, line2), ('Geomagnetic Latitude = 0$^{\circ}$', 'Altitude = 550 km'))
fig.tight_layout()                                # facoltativo: allunga un pelo il grafico
plt.show()
fig.savefig('tesi\Cut-off_latitude_altitude.png', bbox_inches='tight')


figg, ax1g = plt.subplots()
ax1g.set_xlabel('Year')
ax1g.set_ylabel('g$_{1}^{0}$ value')
line1, = ax1g.plot(g, g10(g), color = 'blue')
ax1g.tick_params(axis='x')

ax2g = ax1g.twinx()
line2, = ax2g.plot(g, g10_var(g), color = 'red')
ax2g.set_ylabel('g$_{1}^{0}$ variation')
figg.tight_layout()
ax1g.grid(axis = 'x')
ax1g.grid(axis = 'y')
plt.show()
figg.savefig('tesi\g10 values and variations')

#plt.plot (R, proton_flux(R))

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