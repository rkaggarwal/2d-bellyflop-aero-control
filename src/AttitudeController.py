# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:03:03 2020

@author: raggarwal
"""


import numpy as np
import control
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
from matplotlib import rc
import matplotlib.font_manager
import scipy
from scipy import interpolate
import slycot
from matplotlib import animation
from AeroSensitivityAnalyzer import AeroSensitivityAnalyzer


## ATTITUDE CONTROLLER DESIGN AND SIMULATION

class AttitudeController:
    
    def __init__(self, phi_0 = 0, theta1_0 = 20, theta2_0 = 20, rotInertia = 1):
        
        # Reference position for linearization
        self.phi_0 = phi_0;
        self.theta1_0 = theta1_0;
        self.theta2_0 = theta2_0;
        
        self.J = rotInertia;

        self.aeroAnalyzer = AeroSensitivityAnalyzer();
        
        self.mz_interp_func = self.aeroAnalyzer.getMzInterpFunc();
        self.aeroSensitivities = self.aeroAnalyzer.getAeroSensitivitiesAtConfiguration(self.phi_0, self.theta1_0, self.theta2_0, printOutput = False)
        
        self.dmz_dphi = self.aeroSensitivities[2]*180/np.pi;
        self.dmz_dtheta1 = self.aeroSensitivities[5]*180/np.pi;
        self.dmz_dtheta2 = self.aeroSensitivities[8]*180/np.pi;
        
        
    

    def simulateClosedLoop(self):
        
        n_state = 4;
        n_control = 2;

        theta_offset = np.deg2rad(20); # this is a manual (optional) offset in thetas
        # that trims them to be in a V-shape, which is intuitively more stable in the
        # degrees of freedom not modeled here.  Can be set to 0 and dynamics still work.
        
        q_abs_init = [np.deg2rad(0), np.deg2rad(0), theta_offset, theta_offset];
        
        
        # State: q = [phi; phi_dot; theta1; theta2] (radians, radians/sec)
        # Control: u = [theta1dot; theta2dot] (radians/sec)

        # dynamics is evolved as errors off the reference (which has to be an equilbrium; namely d/dt q_ref = 0)
        # (q - q_ref) = qbar
        # (u - u_ref) = ubar
        # i.e. d/dt (q-q_ref) = A (q-q_ref) + B (u-u_ref)
        # and then at the end the states are offset by the reference values to get
        # the true values.
        
        
        A = np.array([[0, 1, 0, 0],
                      [1/self.J*self.dmz_dphi, 0, 1/self.J*self.dmz_dtheta1, 1/self.J*self.dmz_dtheta2],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0]], dtype = 'float64');
        
        B = np.array([[0, 0],
                      [0, 0],
                      [1, 0],
                      [0, 1]], dtype = 'float64');
        
        # LQR weighting matrices (note the R matrix is a cost on slew rates of 
        # the actuators, while the actual flap angles are in the Q matrix)

        Q = 1*np.diag([30, 15, 1, 1])
        R = 1*np.diag([1, 1])
        
        
       
        
        # solve the infinite-horizon continuous algebraic riccati eqn for X and then solve for K
        (X, L, G) = control.care(A, B, Q, R);
        K = np.linalg.inv(R) @ np.transpose(B) @ X

        sim_time = 20; # seconds
        self.dt = .1; # seconds
        self.time_vector = np.linspace(0, sim_time, int(sim_time/self.dt)+1)
        
        qbar_evolution = np.zeros((n_state, len(self.time_vector))); # radians, radians/sec
        ubar_evolution = np.zeros((n_control, len(self.time_vector))); # radians
        ubar_evolution[:, -1] = [np.nan, np.nan] # control array is one index smaller
        
        self.q_evolution = np.zeros((n_state, len(self.time_vector))); # absolute state
        self.u_evolution = np.zeros((n_control, len(self.time_vector))); # absolute control values
        self.u_evolution[:, -1] = [np.nan, np.nan] # control array is one index smaller
        
        
        self.q_evolution[:, 0] = q_abs_init;
        
        # Set a time-varying reference
        self.q_ref = np.zeros((n_state, len(self.time_vector)))
        
                         
        
        # set body roll references for n distinct periods (can make them all equal if you want)
        roll_refs = np.deg2rad([0, 15, 15, 0, 0, 30, 30, 0, 0]);
        roll_refs_expanded = [ref for ref in roll_refs for i in range(int(1+len(self.time_vector)/len(roll_refs)))]
        self.q_ref[0, :] = roll_refs_expanded[0:len(self.time_vector)]
        
        # body roll rate reference is always 0
        self.q_ref[1, :] = 0;
        
        # solve for u_ref (feedforward)
        for i in range(len(self.time_vector)):
            # solve for the least-norm flap equilibrium angles that can offset aero torque due to the 
            # body roll reference angle.  This amounts to our feed forward command.
            self.q_ref[2, i] = -self.dmz_dphi*self.dmz_dtheta1*self.q_ref[0, i]/(self.dmz_dtheta1**2 + self.dmz_dtheta2**2) + theta_offset;
            self.q_ref[3, i] = -self.dmz_dphi*self.dmz_dtheta2*self.q_ref[0, i]/(self.dmz_dtheta1**2 + self.dmz_dtheta2**2) + theta_offset;
        
        
        # Set the initial conditions
        qbar_evolution[:, 0] = q_abs_init - np.ravel(self.q_ref[:, 0]); 
        
        
        # Simulate closed-loop dynamics
        nonlinearAero = True; # if True, will use real-time interpolation of the CFD nonlinear functions
        
        
        
        for i in range(len(self.time_vector) - 1):
            
            q_curr = self.q_evolution[:, i].reshape(n_state, 1);
            qbar_curr = q_curr - self.q_ref[:, i].reshape(n_state, 1);
            ubar_curr = -K @ qbar_curr
            
            u_curr = ubar_curr; # u_ref == 0 so we don't have to offset anything here
            

            if(nonlinearAero == True):
                delta = 1; # deg, for finite differencing
                
                phi_curr_deg = np.rad2deg(q_curr[0, 0]);
                theta1_curr_deg = np.rad2deg(q_curr[2, 0]);
                theta2_curr_deg = np.rad2deg(q_curr[3, 0]);
                
                dmz_dphi_local = (self.mz_interp_func([phi_curr_deg + delta, theta1_curr_deg, theta2_curr_deg]) - self.mz_interp_func([phi_curr_deg, theta1_curr_deg, theta2_curr_deg]))/delta;
                dmz_dtheta1_local = (self.mz_interp_func([phi_curr_deg, theta1_curr_deg + delta, theta2_curr_deg]) - self.mz_interp_func([phi_curr_deg, theta1_curr_deg, theta2_curr_deg]))/delta;
                dmz_dtheta2_local = (self.mz_interp_func([phi_curr_deg, theta1_curr_deg, theta2_curr_deg + delta]) - self.mz_interp_func([phi_curr_deg, theta1_curr_deg, theta2_curr_deg]))/delta;
                
                avg_dtheta_local = (np.abs(dmz_dtheta1_local) + np.abs(dmz_dtheta2_local))/2;
                dmz_dtheta1_local = np.sign(dmz_dtheta1_local)*avg_dtheta_local;
                dmz_dtheta2_local = np.sign(dmz_dtheta2_local)*avg_dtheta_local;
                
                dmz_dphi_local = dmz_dphi_local[0]*180/np.pi;
                dmz_dtheta1_local = dmz_dtheta1_local[0]*180/np.pi;
                dmz_dtheta2_local = dmz_dtheta2_local[0]*180/np.pi;
                
                
                # update our locally-linear A and B matrices for simulation
                A = np.array([[0, 1, 0, 0],
                          [1/self.J*dmz_dphi_local, 0, 1/self.J*dmz_dtheta1_local, 1/self.J*dmz_dtheta2_local],
                          [0, 0, 0, 0],
                          [0, 0, 0, 0]], dtype = 'float64');
            
                B = np.array([[0, 0],
                          [0, 0],
                          [1, 0],
                          [0, 1]], dtype = 'float64');
        
        
   
            qbar_next = qbar_curr + self.dt*(A @ qbar_curr + B @ ubar_curr)
            
            q_next = qbar_next + self.q_ref[:, i].reshape(n_state, 1);
            
            self.q_evolution[:, i+1] = np.ravel(q_next);
            self.u_evolution[:, i] = np.ravel(u_curr);
            
            qbar_evolution[:, i+1] = np.ravel(qbar_next);
            ubar_evolution[:, i] = np.ravel(ubar_curr);
        
            
        #self.q_evolution[:, -1] = np.ravel(qbar_evolution[:, -1]) + np.ravel(self.q_ref[:, -1])
        
            

        
        fig1, ax1 = plt.subplots(1, 2, figsize=(12, 6));
        
        ax1[0].plot(self.time_vector, np.rad2deg(self.q_evolution[0, :]), label = 'phi');
        ax1[0].plot(self.time_vector, np.rad2deg(self.q_evolution[1, :]), label = 'phi_dot');
        ax1[0].step(self.time_vector, np.rad2deg(self.q_ref[0, :]), 'k--', alpha = .5, label = 'phi ref');
        ax1[0].set_xlabel("Time [s]");
        ax1[0].set_ylabel("Body Roll, Roll Rate [deg, deg/s]");
        ax1[0].set_title("State Variables");
        ax1[0].legend(loc = 'best')
        
        ax1[1].plot(self.time_vector, np.rad2deg(self.q_evolution[2, :]), label = "theta1");
        ax1[1].plot(self.time_vector, np.rad2deg(self.q_evolution[3, :]), label = "theta2");
        ax1[1].set_xlabel("Time [s]");
        ax1[1].set_ylabel("Flap Angle [deg]");
        ax1[1].set_title("Control Variables");
        ax1[1].legend(loc = 'best')
        
        fig1.suptitle("Attitude Reference Tracking, Air Speed = 10 m/s (+y)")
        
        
        
        
        
        
    def animateStateEvolution(self, saveAsFilename = "attitude_ref_tracking_controller"):
    
  
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=10, bitrate = -1)
        
        
        fig2, ax2 = plt.subplots(1, 3, figsize=(18, 6))
        axis_limit = .5;
        ax2[0].set_xlim(-axis_limit, axis_limit)
        ax2[0].set_ylim(-axis_limit, axis_limit)
        ax2[0].set_aspect('equal')
        
        ax2[0].set_xlabel("x-coordinate [m]")
        ax2[0].set_ylabel("y-coordinate [m]");
        ax2[0].set_title("Vehicle State Animation");
        
        
        body_thickness = .05; # m
        body_width= .30; # m
        flap_width = .15; # m


        # Helper function to get the points on the body with attitude changes
        def getBodyPoints(body_roll_angle, theta1_val, theta2_val):
            body_left_point = np.array([[-body_width/2], [0]]);
            body_right_point = np.array([[body_width/2], [0]]);
            
            flap1_tip = body_left_point + np.array([[flap_width*np.cos(np.pi-theta1_val)], [flap_width*np.sin(np.pi-theta1_val)]]);
            flap2_tip = body_right_point + np.array([[flap_width*np.cos(theta2_val)], [flap_width*np.sin(theta2_val)]]);
            
        
            rot_matrix = np.array([[np.cos(body_roll_angle), -np.sin(body_roll_angle)],
                                   [np.sin(body_roll_angle), np.cos(body_roll_angle)]])
        
            
            # note: these rotations only work because the vehicle's center is pinned to the origin
            body_left_point = rot_matrix @ body_left_point;
            body_right_point = rot_matrix @ body_right_point;
            flap1_tip = rot_matrix @ flap1_tip;
            flap2_tip = rot_matrix @ flap2_tip;
            
            return body_left_point, body_right_point, flap1_tip, flap2_tip;
        
        
        
        # Initialize a configuration of patches
        body_line_left = plt.Line2D((0, 0), (-body_width/2, 0), lw = 3, color = 'black') 
        body_line_right = plt.Line2D((0, 0), (body_width/2, 0), lw = 3, color = 'black') 
        
        flap1_line = plt.Line2D((-body_width/2 - flap_width, 0),
                                (-body_width/2, 0), lw = 3, color = 'blue');
             
        flap2_line = plt.Line2D((body_width/2, 0),
                                (body_width/2 + flap_width, 0), lw = 3, color = 'red');                   
             
        
        state_time_line = plt.Line2D((1, 1), (np.rad2deg(np.nanmin(self.q_evolution[0:2, :])), np.rad2deg(np.nanmax(self.q_evolution[0:2, :]))), lw = 2,  color = 'black')
        control_time_line = plt.Line2D((1, 1), (np.rad2deg(np.nanmin(self.q_evolution[2:4, :])), np.rad2deg(np.nanmax(self.q_evolution[2:4, :]))), lw = 2, color = 'black')
        
        
        
        # Plot our state and control time histories for reference
        ax2[1].plot(self.time_vector, np.rad2deg(self.q_evolution[0, :]), label = 'phi');
        ax2[1].plot(self.time_vector, np.rad2deg(self.q_evolution[1, :]), label = 'phi_dot');
        ax2[1].step(self.time_vector, np.rad2deg(self.q_ref[0, :]), 'k--', alpha = .5, label = 'phi ref');
        ax2[1].set_xlabel("Time [s]");
        ax2[1].set_ylabel("Body Roll, Roll Rate [deg, deg/s]");
        ax2[1].set_title("State Variables");
        ax2[1].legend(loc = 'lower right')
        
        ax2[2].plot(self.time_vector, np.rad2deg(self.q_evolution[2, :]), label = "theta1", color = 'blue');
        ax2[2].plot(self.time_vector, np.rad2deg(self.q_evolution[3, :]), label = "theta2", color = 'red');
        ax2[2].set_xlabel("Time [s]");
        ax2[2].set_ylabel("Flap Angle [deg]");
        ax2[2].set_title("Control Variables");
        ax2[2].legend(loc = 'lower right')
        
        fig2.suptitle("Attitude Reference Tracking, Air Speed = 10 m/s (+y)")
            
            
        # Animation functions
        def init():
            ax2[0].add_line(body_line_left)
            ax2[0].add_line(body_line_right)
            ax2[0].add_line(flap1_line)
            ax2[0].add_line(flap2_line)
            
            ax2[1].add_line(state_time_line);
            ax2[2].add_line(control_time_line);

            return []
        
        
        def animationController(i, body_line_left, body_line_right, flap1_line, flap2_line, state_time_line, control_time_line):
            
            animateBodyLineLeft(i, body_line_left)
            animateBodyLineRight(i, body_line_right)
            animateFlap1(i, flap1_line)
            animateFlap2(i, flap2_line)
            
            animateStateTimeLine(i, state_time_line)
            animateControlTimeLine(i, control_time_line)
            
            return [];
        
        
        
        def animateBodyLineLeft(i, line):
            body_roll_angle = self.q_evolution[0, i];
            theta1_val = self.q_evolution[2, i];
            theta2_val = self.q_evolution[3, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([body_left_point[0], 0])
            line.set_ydata([body_left_point[1], 0])
            return line,
        
        
        def animateBodyLineRight(i, line):
            body_roll_angle = self.q_evolution[0, i];
            theta1_val = self.q_evolution[2, i];
            theta2_val = self.q_evolution[3, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([0, body_right_point[0]])
            line.set_ydata([0, body_right_point[1]])
            return line,
        
            
        def animateFlap1(i, line):
            body_roll_angle = self.q_evolution[0, i];
            theta1_val = self.q_evolution[2, i];
            theta2_val = self.q_evolution[3, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([flap1_tip[0], body_left_point[0]])
            line.set_ydata([flap1_tip[1], body_left_point[1]])
            return line,
        
        
        def animateFlap2(i, line):
            body_roll_angle = self.q_evolution[0, i];
            theta1_val = self.q_evolution[2, i];
            theta2_val = self.q_evolution[3, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([body_right_point[0], flap2_tip[0]])
            line.set_ydata([body_right_point[1], flap2_tip[1]])
            return line,
        
        
        def animateStateTimeLine(i, line):
            line.set_xdata([i*self.dt, i*self.dt]);
        
        
        def animateControlTimeLine(i, line):
            line.set_xdata([i*self.dt, i*self.dt]);
        
        anim = animation.FuncAnimation(fig2, animationController, 
                                       init_func=init, 
                                       frames=len(self.time_vector), 
                                       fargs = (body_line_left, body_line_right, flap1_line, flap2_line, state_time_line, control_time_line),
                                       interval=self.dt*1000,
                                       blit=True)
        
        anim.save(saveAsFilename + '.gif', writer=writer)
        
        plt.show()    
        
    