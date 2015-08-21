#!/usr/bin/env python

import sys
import randomfield

generator = randomfield.Generator(8, 256, 4096, grid_spacing_Mpc_h=1.0, verbose=True)
delta = generator.generate_delta_field(smoothing_length_Mpc_h=2.0, seed=123, show_plot=True)
