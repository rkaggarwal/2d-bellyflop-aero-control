# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:02:20 2020

@author: raggarwal
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate


## CFD POST-PROCESSING AND ANALYSIS

class AeroSensitivityAnalyzer:
   
    def __init__(self, filename = "cfd_parameter_study/full_cfd_sensitivity_runs.csv"):
        
        self.filename = filename;
        
        # CSV column ordering: Run Name, air speed (m/s), phi_from_+y (deg), theta1_from_-x (deg, clockwise is positive),
        # theta2_from_+x (deg, counter clockwise is positive), fx_net (N), fy_net (N), mz_net (N-m)
        self.cfdData = pd.read_csv(self.filename)


        
        # Generate interpolating functions from discrete CFD data        
        # Note, these lists are hard-coded based off of what was run in CFD
        self.phi_list = np.array([-45, -30, -15, 0, 15, 30, 45], dtype = 'float64'); # deg
        self.theta1_list = np.array([-60, -40, -20, 0,  20, 40, 60], dtype = 'float64'); # deg
        self.theta2_list = np.array([-60, -40, -20, 0, 20, 40, 60], dtype = 'float64'); # deg
        

        # for the following grids:
        # the 1st coord of the grid is the phi_val, 2nd coord is theta1_val, 3rd coord is theta2_val
        self.phi_grid, self.theta1_grid, self.theta2_grid = np.meshgrid(self.phi_list, self.theta1_list, self.theta2_list, indexing = 'ij');
        
        
        self.fx_grid = np.zeros((len(self.phi_list), len(self.theta1_list), len(self.theta2_list)));
        self.fy_grid = np.zeros((len(self.phi_list), len(self.theta1_list), len(self.theta2_list)));
        self.mz_grid = np.zeros((len(self.phi_list), len(self.theta1_list), len(self.theta2_list)));
        
        
        
        for i in range(len(self.phi_list)):
            phi_val = self.phi_list[i];
            phi_val_cfd = self.phi_list[i]; # this is the phi val from the raw .csv and is the csys in CFD  
            if(phi_val_cfd == 0): phi_val_cfd = 0.1; # the .1 values were a hack to get the DesignModeler CAD to work for 0 deg
        
            
            for j in range(len(self.theta1_list)):
                theta1_val = self.theta1_list[j];
                theta1_val_cfd = 90 - theta1_val;
                
                for k in range(len(self.theta2_list)):
                    theta2_val = self.theta2_list[k];
                    theta2_val_cfd = 90 - theta2_val;
                    
                    
                    
                    if(phi_val >= 0): # this means we have direct data for this
                        temp = self.cfdData.loc[(self.cfdData['phi_from_y_deg'] == phi_val_cfd) &
                                            (self.cfdData['theta1_from_y_deg'] == theta1_val_cfd) &
                                            (self.cfdData['theta2_from_y_deg'] == theta2_val_cfd)]
        
                        self.fx_grid[i, j, k] = temp.iloc[0, 5];
                        self.fy_grid[i, j, k] = temp.iloc[0, 6];
                        self.mz_grid[i, j, k] = temp.iloc[0, 7];
                    
                    else: # we have to leverage symmetry in the solutions (this symmetry was leveraged to reduce CFD run times)
                        temp = self.cfdData.loc[(self.cfdData['phi_from_y_deg'] == -1*phi_val_cfd) &
                                            (self.cfdData['theta1_from_y_deg'] == theta2_val_cfd) &
                                            (self.cfdData['theta2_from_y_deg'] == theta1_val_cfd)]
        
                        # fx and mz have opposite signs, fy stays the same (drag)
                        self.fx_grid[i, j, k] = -temp.iloc[0, 5];
                        self.fy_grid[i, j, k] = temp.iloc[0, 6];
                        self.mz_grid[i, j, k] = -temp.iloc[0, 7];
        
        self.fx_interp_func = interpolate.RegularGridInterpolator((self.phi_list, self.theta1_list, self.theta2_list), self.fx_grid)
        self.fy_interp_func = interpolate.RegularGridInterpolator((self.phi_list, self.theta1_list, self.theta2_list), self.fy_grid)
        self.mz_interp_func = interpolate.RegularGridInterpolator((self.phi_list, self.theta1_list, self.theta2_list), self.mz_grid)
        
        
        
    def getFxInterpFunc(self):
        return self.fx_interp_func;
    
    def getFyInterpFunc(self):
        return self.fy_interp_func;
    
    def getMzInterpFunc(self):
        return self.mz_interp_func;
    
    
    
    # phi, theta1, theta2 have to all be in degrees
    def getAeroSensitivitiesAtConfiguration(self, phi_0 = 0, theta1_0 = 20, theta2_0 = 20, printOutput = True):
    

        delta = 1; # deg, for finite differencing

    
        dfx_dphi = (self.fx_interp_func([phi_0 + delta, theta1_0, theta2_0]) - self.fx_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dfy_dphi = (self.fy_interp_func([phi_0 + delta, theta1_0, theta2_0]) - self.fy_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dmz_dphi = (self.mz_interp_func([phi_0 + delta, theta1_0, theta2_0]) - self.mz_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        
        dfx_dtheta1 = (self.fx_interp_func([phi_0, theta1_0 + delta, theta2_0]) - self.fx_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dfy_dtheta1 = (self.fy_interp_func([phi_0, theta1_0 + delta, theta2_0]) - self.fy_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dmz_dtheta1 = (self.mz_interp_func([phi_0, theta1_0 + delta, theta2_0]) - self.mz_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        
        dfx_dtheta2 = (self.fx_interp_func([phi_0, theta1_0, theta2_0 + delta]) - self.fx_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dfy_dtheta2 = (self.fy_interp_func([phi_0, theta1_0, theta2_0 + delta]) - self.fy_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        dmz_dtheta2 = (self.mz_interp_func([phi_0, theta1_0, theta2_0 + delta]) - self.mz_interp_func([phi_0, theta1_0, theta2_0]))/delta;
        
        if(printOutput):
            print("At vehicle state (phi, theta1, theta2) = ({}, {}, {}) deg,\n".format(phi_0, theta1_0, theta2_0),
                  "dfx_dphi: {} N/deg \n".format(dfx_dphi),
                  "dfy_dphi: {} N/deg \n".format(dfy_dphi),
                  "dmz_dphi: {} N-m/deg \n\n".format(dmz_dphi),
                  "dfx_dtheta1: {} N/deg \n".format(dfx_dtheta1),
                  "dfy_dtheta1: {} N/deg \n".format(dfy_dtheta1),
                  "dmz_dtheta1: {} N-m/deg \n\n".format(dmz_dtheta1),
                  "dfx_dtheta2: {} N/deg \n".format(dfx_dtheta2),
                  "dfy_dtheta2: {} N/deg \n".format(dfy_dtheta2),
                  "dmz_dtheta2: {} N-m/deg \n".format(dmz_dtheta2))
    
    
        return [dfx_dphi, dfy_dphi, dmz_dphi,
                dfx_dtheta1, dfy_dtheta1, dmz_dtheta1,
                dfx_dtheta2, dfy_dtheta2, dmz_dtheta2]
    
    
    
    def plotZeroMomentPoses(self, moment_threshold = .1):
        # Find all points that have near-zero moment and scatter plot them
        # phi, theta1 and theta2
        

        soln_points = np.zeros((3, len(self.phi_list)*len(self.theta1_list)*len(self.theta2_list)))
        soln_counter = 0;
        
        
        for i in range(len(self.phi_list)):
            phi_val = self.phi_list[i];
            
            for j in range(len(self.theta1_list)):
                theta1_val = self.theta1_list[j];
                
                for k in range(len(self.theta2_list)):
                    theta2_val = self.theta2_list[k];
                    
                    if(np.abs(self.mz_grid[i, j, k]) < moment_threshold):
                        
                        soln_points[:, soln_counter] = [phi_val, theta1_val, theta2_val];
                        soln_counter += 1;
                        
        
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_trisurf(soln_points[0, :], soln_points[1, :], soln_points[2, :], color = 'b', edgecolor = 'k', alpha = .5)
        ax.set_xlabel("Body Roll (phi) [deg]");
        ax.set_ylabel("Flap 1 Angle (theta2) [deg]")
        ax.set_zlabel("Flap 2 Angle (theta2) [deg]")
        ax.set_title("Zero Moment Poses")
        plt.tight_layout()
        return soln_points;
    
    
    def plotFullSensitivities(self):

        
        fig0, ax0 = plt.subplots(3, 3, constrained_layout=True);
        num_lines = 75; # number of contour lines in an individual plot
        mdpt = int(len(self.phi_list)/2)
        
        
        
        # First row is Fx sensitivities
        fx_min = np.min(self.fx_grid);
        fx_max = np.max(self.fx_grid);
        
        cb00 = ax0[0, 0].contourf(self.theta1_grid[0, :, :], self.theta2_grid[0, :, :], self.fx_grid[0, :, :], levels = np.linspace(fx_min, fx_max, num_lines), cmap = 'jet')
        ax0[0, 0].set_xlabel('theta1 [deg]')
        ax0[0, 0].set_ylabel('theta2 [deg]')
        ax0[0, 0].set_title("Fx vs. Flap Angles, phi = {} deg".format(self.phi_list[0]))
        cbar00 = fig0.colorbar(cb00, ax = ax0[0, 0])
        cbar00.ax.set_ylabel('Fx [N]')
        
        cb01 = ax0[0, 1].contourf(self.theta1_grid[mdpt, :, :], self.theta2_grid[mdpt, :, :], self.fx_grid[mdpt, :, :], levels = np.linspace(fx_min, fx_max, num_lines), cmap = 'jet')
        ax0[0, 1].set_xlabel('theta1 [deg]')
        ax0[0, 1].set_ylabel('theta2 [deg]')
        ax0[0, 1].set_title("Fx vs. Flap Angles, phi = {} deg".format(self.phi_list[mdpt]))
        cbar01 = fig0.colorbar(cb01, ax = ax0[0, 1])
        cbar01.ax.set_ylabel('Fx [N]')
        cb01.set_clim(vmin = np.min(self.fx_grid), vmax = np.max(self.fx_grid))
        
        cb02 = ax0[0, 2].contourf(self.theta1_grid[-1, :, :], self.theta2_grid[-1, :, :], self.fx_grid[-1, :, :], levels = np.linspace(fx_min, fx_max, num_lines), cmap = 'jet')
        ax0[0, 2].set_xlabel('theta1 [deg]')
        ax0[0, 2].set_ylabel('theta2 [deg]')
        ax0[0, 2].set_title("Fx vs. Flap Angles, phi = {} deg".format(self.phi_list[-1]))
        cbar02 = fig0.colorbar(cb02, ax = ax0[0, 2])
        cbar02.ax.set_ylabel('Fx [N]')
        
        
        
        
        # Second row is Fy sensitivities
        fy_min = np.min(self.fy_grid);
        fy_max = np.max(self.fy_grid);
        
        cb10 = ax0[1, 0].contourf(self.theta1_grid[0, :, :], self.theta2_grid[0, :, :], self.fy_grid[0, :, :], levels = np.linspace(fy_min, fy_max, num_lines), cmap = 'jet')
        ax0[1, 0].set_xlabel('theta1 [deg]')
        ax0[1, 0].set_ylabel('theta2 [deg]')
        ax0[1, 0].set_title("Fy vs. Flap Angles, phi = {} deg".format(self.phi_list[0]))
        cbar10 = fig0.colorbar(cb10, ax = ax0[1, 0])
        cbar10.ax.set_ylabel('Fy [N]')
        
        cb11 = ax0[1, 1].contourf(self.theta1_grid[mdpt, :, :], self.theta2_grid[mdpt, :, :], self.fy_grid[mdpt, :, :], levels = np.linspace(fy_min, fy_max, num_lines), cmap = 'jet')
        ax0[1, 1].set_xlabel('theta1 [deg]')
        ax0[1, 1].set_ylabel('theta2 [deg]')
        ax0[1, 1].set_title("Fy vs. Flap Angles, phi = {} deg".format(self.phi_list[mdpt]))
        cbar11 = fig0.colorbar(cb11, ax = ax0[1, 1])
        cbar11.ax.set_ylabel('Fy [N]')
        
        cb12 = ax0[1, 2].contourf(self.theta1_grid[-1, :, :], self.theta2_grid[-1, :, :], self.fy_grid[-1, :, :], levels = np.linspace(fy_min, fy_max, num_lines), cmap = 'jet')
        ax0[1, 2].set_xlabel('theta1 [deg]')
        ax0[1, 2].set_ylabel('theta2 [deg]')
        ax0[1, 2].set_title("Fy vs. Flap Angles, phi = {} deg".format(self.phi_list[-1]))
        cbar12 = fig0.colorbar(cb12, ax = ax0[1, 2])
        cbar12.ax.set_ylabel('Fy [N]')
        
        
        
        
        # Third row is the Mz sensitivities
        mz_min = np.min(self.mz_grid);
        mz_max = np.max(self.mz_grid);
        
        cb20 = ax0[2, 0].contourf(self.theta1_grid[0, :, :], self.theta2_grid[0, :, :], self.mz_grid[0, :, :], levels = np.linspace(mz_min, mz_max, num_lines), cmap = 'jet')
        ax0[2, 0].set_xlabel('theta1 [deg]')
        ax0[2, 0].set_ylabel('theta2 [deg]')
        ax0[2, 0].set_title("Mz vs. Flap Angles, phi = {} deg".format(self.phi_list[0]))
        cbar20 = fig0.colorbar(cb20, ax = ax0[2, 0])
        cbar20.ax.set_ylabel('Mz [N-m]')
        
        cb21 = ax0[2, 1].contourf(self.theta1_grid[mdpt, :, :], self.theta2_grid[mdpt, :, :], self.mz_grid[mdpt, :, :], levels = np.linspace(mz_min, mz_max, num_lines), cmap = 'jet')
        ax0[2, 1].set_xlabel('theta1 [deg]')
        ax0[2, 1].set_ylabel('theta2 [deg]')
        ax0[2, 1].set_title("Mz vs. Flap Angles, phi = {} deg".format(self.phi_list[mdpt]))
        cbar21 = fig0.colorbar(cb21, ax = ax0[2, 1])
        cbar21.ax.set_ylabel('Mz [N-m]')
        
        cb22 = ax0[2, 2].contourf(self.theta1_grid[-1, :, :], self.theta2_grid[-1, :, :], self.mz_grid[-1, :, :], levels = np.linspace(mz_min, mz_max, num_lines), cmap = 'jet')
        ax0[2, 2].set_xlabel('theta1 [deg]')
        ax0[2, 2].set_ylabel('theta2 [deg]')
        ax0[2, 2].set_title("Mz vs. Flap Angles, phi = {} deg".format(self.phi_list[-1]))
        cbar22 = fig0.colorbar(cb22, ax = ax0[2, 2])
        cbar22.ax.set_ylabel('Mz [N-m]')
        
        fig0.suptitle("Aerodynamic Load Sensitivities, u = 10 m/s")



