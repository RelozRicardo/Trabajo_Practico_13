#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp
from sympy import symbols, Matrix
from sympy import init_printing
from sympy import simplify
from sympy import symbols, expand, factor, solve
from pytc2.general import to_latex

from schemdraw import Drawing

# Ahora importamos las funciones de PyTC2

from pytc2.remociones import remover_polo_dc, remover_polo_jw , remover_polo_jw2 ,remover_polo_infinito , remover_polo_dc2 , remover_polo_infinito2, isFRP , remover_polo_sigma
from pytc2.dibujar import display, dibujar_tanque_serie, dibujar_puerto_entrada, dibujar_funcion_exc_abajo,  dibujar_elemento_serie, dibujar_elemento_derivacion,  dibujar_tanque_derivacion, dibujar_tanque_RC_serie,  dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC
from pytc2.dibujar import dibujar_Pi, dibujar_Tee, dibujar_lattice
from pytc2.sintesis_dipolo import cauer_LC

from pytc2.dibujar import dibujar_cauer_LC
from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s, simplify_n_monic
from IPython.display import display,  Markdown

from pytc2.cuadripolos import Z2Tabcd_s, Y2Tabcd_s, Tabcd2Z_s, Tabcd2Y_s
from pytc2.cuadripolos import calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn

import scipy.signal as sig
from pytc2.sistemas_lineales import analyze_sys
from sympy import conjugate, I,pi, N

patron = "\n" + "/" * 75 + "\n" # 75 barras para ajustar el largo deseado
# Activar la impresión en formato LaTeX
init_printing()
# Definir la variable simbólica s
s = sp.symbols('s',imaginary=True)
w = sp.symbols('w', real = True)
# Definir la función de transferencia T(s)
numerador = 15
denominador = s**3 + 6*s**2 + 15*s + 15
Tb = numerador / denominador
Tbc = sp.conjugate(Tb)

print_latex(a_equal_b_latex_s('Tb(s)', Tb))
print_latex(a_equal_b_latex_s('Tbc(s)', Tbc))

print(patron)

S21 = Tb
S21m2 = S21 * sp.conjugate(S21)
S21m2 = sp.simplify(sp.expand(sp.factor(S21m2)))
print_latex(a_equal_b_latex_s('S21m2(s)', S21m2))

print(patron)

S11m2 = 1 - S21m2
print_latex(a_equal_b_latex_s('S11m2(s)', S11m2))


S11m2 = sp.simplify(sp.expand(sp.factor(S11m2)))
print_latex(a_equal_b_latex_s('S11m2(s)', S11m2))

#//////////////////////////////////////////////////////////////////////////////
#s = sp.symbols('s',Complex=True)

print(patron)
A = sp.sqrt(6 + 2*sp.sqrt(45)) ; B = sp.sqrt(45) ; 
D = 6 ; E = 15 ; F = 15

NUMs11 =s**3 + A*s**2 + B*s 
DENs11 = s**3 + D*s**2 + E*s + F
S11 = NUMs11/DENs11
print_latex(a_equal_b_latex_s('S11(s)', S11))
print_latex(a_equal_b_latex_s('S11(s)', S11.evalf(2)))


print("\nVERIFICACION DE S11 obteniendo el Modulo cuadrado")
print_latex(a_equal_b_latex_s('S11m2(s)', sp.simplify(sp.expand(sp.factor(S11 * sp.conjugate(S11))))))

Z1 = (1+S11)/(1-S11)
Z1 = sp.simplify(sp.expand(sp.factor(Z1)))


H = sig.TransferFunction( [2, 10.4064053224, 21.7082039325 , 15], [1.59359467763, 8.2917960675, +15] )

analyze_sys(H)

print_latex(a_equal_b_latex_s('Z1(s)', Z1))
#//////////////////////////////////////////////////////////////////////////////
print(patron)
print("Modificacion de circuito, con R = 50 ohm y wo = 2*pi*10**6")

print(patron)
##/////////////////////////////////////////////////////////////////////////////

wo = 2*sp.pi*10**6
rg = 50

numerador = 15
denominador = (s/wo)**3 + 6*(s/wo)**2 + 15*(s/wo) + 15

S11 = sp.simplify(S11.subs(s, s/wo))
print_latex(a_equal_b_latex_s('S11(s)', S11))
print_latex(a_equal_b_latex_s('S11(s)', S11.evalf(2)))

Z1 = rg*(1+S11)/(1-S11)
Z1 = sp.simplify(sp.expand(sp.factor(Z1)))
#Z1 = sp.nsimplify(sp.expand(sp.factor(Z1))) #(para saber el desarme)
print_latex(a_equal_b_latex_s('Z1(s)', Z1))
print_latex(a_equal_b_latex_s('Z1(s)', Z1.evalf(2)))

