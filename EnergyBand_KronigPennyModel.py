# Energy Bands of the Kronig Penny Model graphical representation
# Author: Jonas Hoecker

import numpy as np
import math as m
import matplotlib.pyplot as plt

v = 10
def equation(e, cos_kxi):
    return m.cos(m.sqrt(e)) * m.cos(m.sqrt(v + e)) - ((v/2) + e)/(m.sqrt(e * (v + e))) * m.sin(m.sqrt(e)) * m.sin(m.sqrt(v + e)) - cos_kxi

def numerical_derivative(f, x, cos_kxi, h=1e-6):
    return (f(x + h, cos_kxi) - f(x, cos_kxi))/h

def newton(f, initial_guess, cos_kxi, tolerance=1e-6):
    x = initial_guess
    ic = 0
    d = 0.
    der_f = 0.
    fn = 0.
    while True:
        ic += 1
        fn = f(x, cos_kxi)
        der_f = numerical_derivative(f, x, cos_kxi)
        new_x = x - fn / der_f
        d = m.fabs(new_x - x)
        print (x,fn,der_f)
        if m.fabs(der_f) < 0.0001 or ic > 100:
            print ('we are here',der_f,ic)
            x += 0.1
            break
        if m.fabs(new_x - x) < tolerance:
            print('got it')
            x = new_x
            break
        x = new_x
    print (x,d,fn,der_f,ic,cos_kxi)
    return x


def find_e_values(guess, equation):
    e_values = []
    kxi_values = []
    for kxi in np.arange(0., 2*m.pi, 0.001):
        cos_kxi = m.cos(kxi)
        print('-------------------------------', kxi, cos_kxi)
        e_value = newton(equation, guess, cos_kxi)
        e_values.append(e_value)
        kxi_values.append(kxi)
        guess = e_value
    return kxi_values, e_values

guess_values = [6, 34, 84]
for index, guess in enumerate(guess_values):
    kxi_values, e_values = find_e_values(guess, equation)
    plt.plot(kxi_values, e_values, label=f"Energy band {index+1}")

plt.title("Energy bands of the Kronig Penney Model")
plt.xlabel(r'$k\xi$')
plt.ylabel(r'Dimensionless energy $\varepsilon$')
plt.grid(True)
plt.legend()
plt.show()