# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:46:56 2020

@author: raggarwal
"""

from AttitudeController import AttitudeController
from AeroSensitivityAnalyzer import AeroSensitivityAnalyzer;
from HorizonOptimization import HorizonOptimization;
import matplotlib.pyplot as plt


plt.close('all')

# Plot aero sensitivities from CFD
# aeroSens = AeroSensitivityAnalyzer();
# aeroSens.plotFullSensitivities();

# Simulate closed-loop attitude reference tracking
at_cont = AttitudeController(phi_0 = 0, theta1_0 = 20, theta2_0 = 20, rotInertia = 1);
at_cont.simulateClosedLoop();
# at_cont.animateStateEvolution(saveAsFilename = "attitude_ref_tracking_video");

# Plan a non-linear single horizon
# horzOpt = HorizonOptimization(phi_0 = 0, theta1_0 = 20, theta2_0 = 20, rotInertia = 1, mass = 1);
# horzOpt.planHorizon(x_init = 100, y_init = 100, plottingOn = True);
# horzOpt.animateHorizon(saveAsFilename = "minimal_cost_trajectory_video");
 