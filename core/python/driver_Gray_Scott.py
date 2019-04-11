#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 09:00:21 2017

@author: blake

Gray Scott system.

"""


# -----------------------------------------------------------------------------
# imports
# -----------------------------------------------------------------------------

from sympy import symbols
from time_evolution import create_code

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

# parameters
alpha, gamma = symbols('alpha, gamma')

# system variables
u,v = symbols('u,v')

# vector containing the system variables
U = [u,v]

# system parameters
parameters = [alpha, gamma]

#
# f_0(U)_t + F_1(U)_x + G(U) = (B(U)U_x)_x 
#

f0 = [u,v]

f1 = [0,0]

BU = [[1,0],[0,1]]

G = [(u*v**2-alpha*(1-u)),(-u*v**2/gamma +v/gamma)]


# file path to location where code will be written
file_path = '/Users/blake/Dropbox/stablab20/gray_scott/gray_scott_code'

create_code(U,parameters,f0,f1,G,BU,file_path)

