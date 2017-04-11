#!/usr/bin/env python
"""Wrapper of WenOOF Fortran library"""

from __future__ import print_function
from ctypes import CDLL, c_char_p, c_double, c_int, POINTER
import matplotlib.pyplot as plt
import numpy as np

stencils = 3
interpolator_type = 'interpolator-JS\n'
eps = 1e-19
x_target = 0.3
verbose = 1

wenoof = CDLL('./shared/libwenoof.so')
wenoof.wenoof_initialize_c_wrap.argtypes = [c_char_p, c_int, c_double, c_double, c_int]
wenoof.wenoof_interpolate_c_wrap.argtypes = [c_int, POINTER(c_double), POINTER(c_double)]
wenoof.wenoof_reconstruct_c_wrap.argtypes = [c_int, POINTER(c_double), POINTER(c_double)]

wenoof.wenoof_initialize_c_wrap(interpolator_type, stencils, eps, x_target, verbose)

cells_number = 50
x_cell, dx = np.linspace(start=-np.pi, stop=np.pi, num=cells_number + 2 * stencils, endpoint=True, retstep=True)
y_cell = np.sin(x_cell)
y_weno = y_cell
for i, x in enumerate(x_cell):
  if i >= stencils and i < cells_number + stencils:
    print(i, 'in')
    wenoof.wenoof_interpolate_c_wrap(stencils,
                                     y_cell[i-stencils:].ctypes.data_as(POINTER(c_double)),
                                     y_weno[i:i].ctypes.data_as(POINTER(c_double)))
  else:
    print(i, 'out')

x_interp = x_cell + x_target * dx
y_interp = np.sin(x_interp)
plt.plot(x_cell[stencils:], y_cell[stencils:])
plt.plot(x_cell[stencils:], y_weno[stencils:], 'ro', ms=1.5)
plt.xlabel('Angle [rad]')
plt.ylabel('sin(x)')
plt.axis('tight')
plt.show()
