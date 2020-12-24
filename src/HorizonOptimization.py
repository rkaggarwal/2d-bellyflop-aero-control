# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 22:20:18 2020

@author: raggarwal
"""

from casadi import *
import matplotlib.pyplot as plt
import numpy as np
from AeroSensitivityAnalyzer import AeroSensitivityAnalyzer;
from matplotlib import animation;
    


## TRAJECTORY OPTIMIZATION


class HorizonOptimization:
    
    def __init__(self, phi_0 = 0, theta1_0 = 20, theta2_0 = 20, rotInertia = 1, mass = 1):
        
        # Reference position for linearization
        self.phi_0 = phi_0;
        self.theta1_0 = theta1_0;
        self.theta2_0 = theta2_0;
        
        self.J = rotInertia;
        self.mass = mass;
        
        self.aeroAnalyzer = AeroSensitivityAnalyzer();
        
        
        self.fx_interp_func = self.aeroAnalyzer.getFxInterpFunc();
        self.fy_interp_func = self.aeroAnalyzer.getFyInterpFunc();
        self.mz_interp_func = self.aeroAnalyzer.getMzInterpFunc();
        self.aeroSensitivities = self.aeroAnalyzer.getAeroSensitivitiesAtConfiguration(self.phi_0, self.theta1_0, self.theta2_0, printOutput = False)
        
        
        # loads at the linearization point
        self.fx_const = 0;
        self.fy_const = self.fy_interp_func([self.phi_0, self.theta1_0, self.theta2_0])[0];
        self.mz_const = 0;
        
        # sensitivities are returned in N/deg or N-m/deg and converted to N/rad or N-m/rad
        self.dfx_dphi = self.aeroSensitivities[0]*180/np.pi;
        self.dfx_dtheta1 = self.aeroSensitivities[3]*180/np.pi;
        self.dfx_dtheta2 = self.aeroSensitivities[6]*180/np.pi;
        
        self.dfy_dphi = self.aeroSensitivities[1]*180/np.pi;
        self.dfy_dtheta1 = self.aeroSensitivities[4]*180/np.pi;
        self.dfy_dtheta2 = self.aeroSensitivities[7]*180/np.pi;
        
        
        self.dmz_dphi = self.aeroSensitivities[2]*180/np.pi;
        self.dmz_dtheta1 = self.aeroSensitivities[5]*180/np.pi;
        self.dmz_dtheta2 = self.aeroSensitivities[8]*180/np.pi;
        
        
        self.n_states = 9; # X, Y, X_dot, Y_dot, phi, phi_dot, theta1, theta2, g
        self.n_control = 2; # theta1_dot, theta2_dot



        # Linearized A and B matrices, y accels can be nullified out
        # self.A = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
        #                    [0, 0, 0, 0, self.dfx_dphi/self.mass, 0, self.dfx_dtheta1/self.mass, self.dfx_dtheta2/self.mass, 0],
        #                    [0, 0, 0, 1, 0, 0, 0, 0, 0],
        #                    [0, 0, 0, 0, 0, 0, 0, 0, 0], #[0, 0, 0, 0, self.dfy_dphi/self.mass, 0, self.dfy_dtheta1/self.mass, self.dfy_dtheta2/self.mass, -1],
        #                    [0, 0, 0, 0, 0, 1, 0, 0, 0],
        #                    [0, 0, 0, 0, self.dmz_dphi/self.J, 0, self.dmz_dtheta1/self.J, self.dmz_dtheta2/self.J, 0],
        #                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
        #                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
        #                    [0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype = 'float');
        
        # self.B = np.array([[0, 0],
        #                    [0, 0],
        #                    [0, 0],
        #                    [0, 0],
        #                    [0, 0],
        #                    [0, 0],
        #                    [1, 0],
        #                    [0, 1],
        #                    [0, 0]], dtype = 'float')
                           






    def dq_dt(self, q, u):
        """
        Non-linear state derivative
        """
        
        # Aero coefficients are scaled by the y velocity according to a
        # standard square law for aerodynamic loads
        
        K_y_vel = (q[3]/10)**2 # aero coeff's generated at 10 m/s from CFD
        
        x = q[0];
        xdot = q[1];
        y = q[2];
        ydot = q[3];
        phi = q[4];
        phidot = q[5];
        theta1 = q[6];
        theta2 = q[7];
        g = q[8];
        
        # Non-linear state derivatives
        x_dot = xdot;
        xdot_dot = K_y_vel*self.dfx_dphi/self.mass*(phi-np.deg2rad(self.phi_0)) + K_y_vel*self.dfx_dtheta1/self.mass*(theta1-np.deg2rad(self.theta1_0)) + K_y_vel*self.dfx_dtheta2/self.mass*(theta2-np.deg2rad(self.theta2_0));
        y_dot = ydot;
        ydot_dot = K_y_vel*self.fy_const/self.mass + K_y_vel*self.dfy_dphi/self.mass*(phi-np.deg2rad(self.phi_0)) + K_y_vel*self.dfy_dtheta1/self.mass*(theta1-np.deg2rad(self.theta1_0)) + K_y_vel*self.dfy_dtheta2/self.mass*(theta2-np.deg2rad(self.theta2_0)) - g/self.mass;
        phi_dot = phidot;
        phidot_dot = K_y_vel*self.dmz_dphi/self.J*(phi-np.deg2rad(self.phi_0)) + K_y_vel*self.dmz_dtheta1/self.J*(theta1-np.deg2rad(self.theta1_0)) + K_y_vel*self.dmz_dtheta2/self.J*(theta2-np.deg2rad(self.theta2_0));
        theta1_dot = u[0];
        theta2_dot = u[1];
        g_dot = 0;
        
        return vertcat(x_dot, xdot_dot, y_dot, ydot_dot, phi_dot, phidot_dot, theta1_dot, theta2_dot, g_dot)


    
    
    
    
    def planHorizon(self, x_init = 50, y_init = 50, plottingOn = True):
        """
        Plan a single-shot trajectory.  Note this is a bit crude,
        there's really a lot more out there on how to do this
        more efficiently/in a principled way
        """
        
        N = 50; # number of discretization points
        
        opti = Opti();
        
        Q = opti.variable(self.n_states, N+1); # state trajectory
        
        X = Q[0, :];
        X_DOT = Q[1, :];
        Y = Q[2, :];
        Y_DOT = Q[3, :];
        PHI = Q[4, :];
        PHI_DOT = Q[5, :]
        THETA1 = Q[6, :];
        THETA2 = Q[7, :];
        G = Q[8, :];
        
        U = opti.variable(self.n_control, N); # control trajectory
        THETA1_DOT = U[0, :];
        THETA2_DOT = U[1, :];
        
        T = opti.variable(); # s, final time
        
        DT = T/N; # s, time interval
        
        #opti.minimize(T**2); # min time trajectories
        #opti.minimize(Y_DOT[-1]**2) # min final y velocity trajectories
        
        
        PHI_DOT_NORM = 0;
        THETA_DOT_NORM = 0;
        for ii in range(N):
            PHI_DOT_NORM += PHI_DOT[ii]**2;
            THETA_DOT_NORM += THETA1_DOT[ii]**2 + THETA2_DOT[ii]**2;
        
        # min body roll, final time and actuator effort trajectories
        # (can set different relative weights for experimentation)
        opti.minimize(1*PHI_DOT_NORM + 1*T**2 + 1*THETA_DOT_NORM)
        
        
        for ii in range(N):
            k1 = self.dq_dt(Q[:, ii],              U[:, ii]);
            k2 = self.dq_dt(Q[:, ii] + DT/2*k1,    U[:, ii]);
            k3 = self.dq_dt(Q[:, ii] + DT/2*k2,    U[:, ii]);
            k4 = self.dq_dt(Q[:, ii] + DT*k3,      U[:, ii]);
            
            q_next = Q[:, ii] + DT/6*(k1 + 2*k2 + 2*k3 + k4);
        
            opti.subject_to(Q[:, ii+1] == q_next); # RK4 integration
            #opti.subject_to(Q[:, ii+1] == Q[:, ii] + DT*dq_dt(Q[:, ii], U[:, ii])); # Euler integration
        
        
        # vehicle configuration constraints
        opti.subject_to(opti.bounded(-np.deg2rad(45), PHI, np.deg2rad(45)));
        opti.subject_to(opti.bounded(-np.deg2rad(45), THETA1, np.deg2rad(45)));
        opti.subject_to(opti.bounded(-np.deg2rad(45), THETA2, np.deg2rad(45)));
        
        # flap slew rate constraints
        opti.subject_to(opti.bounded(-np.deg2rad(30), THETA1_DOT, np.deg2rad(30)));
        opti.subject_to(opti.bounded(-np.deg2rad(30), THETA2_DOT, np.deg2rad(30)));
        
        # misc constraints
        opti.subject_to(T > 0);
        
        
        # initial conditions
        opti.subject_to(X[0] == x_init);
        opti.subject_to(X_DOT[0] == 0);
        opti.subject_to(Y[0] == y_init);
        opti.subject_to(Y_DOT[0] == -10); # m/s
        opti.subject_to(PHI[0] == 0);
        opti.subject_to(PHI_DOT[0] == 0);
        opti.subject_to(THETA1[0] == np.deg2rad(20));
        opti.subject_to(THETA2[0] == np.deg2rad(20));
        opti.subject_to(G == 9.81); # m/s^2
        
        # terminal conditions
        opti.subject_to(X[-1] == 0);
        opti.subject_to(X_DOT[-1] == 0);
        opti.subject_to(Y[-1] == 0);
        # constraint on final y velocity
        #opti.subject_to(Y_DOT[-1] == 0);
        opti.subject_to(PHI[-1] == 0);
        opti.subject_to(PHI_DOT[-1] == 0);
        #opti.subject_to(THETA1[-1] == np.deg2rad(20));
        #opti.subject_to(THETA2[-1] == np.deg2rad(20));
        #opti.subject_to(THETA1[-1] > 0);
        #opti.subject_to(THETA2[-1] > 0);
        opti.subject_to(THETA1[-1] == THETA2[-1])

        opti.set_initial(T, 5)
        
        
        
        p_opts = {"expand": True}
        s_opts = {"print_level" : 5,  "acceptable_constr_viol_tol": 1e-3, "constr_viol_tol": 1e-3, "acceptable_tol" : 1e-3, "acceptable_obj_change_tol": 1e-3}
        opti.solver("ipopt", p_opts, s_opts)
        sol = opti.solve()
        
        Q_sol = sol.value(Q);
        U_sol = sol.value(U);
        T_sol = sol.value(T);
        DT_sol = sol.value(DT);
        T_VECTOR_sol = np.linspace(0, T_sol, N+1)
        
        
        self.Q_sol = Q_sol;
        self.T_VECTOR_sol = T_VECTOR_sol;
        self.DT = DT_sol;
        
        
        
        if(plottingOn):
            fig0, ax0 = plt.subplots(1, 3, figsize=(18, 6));
            fig0.suptitle("Minimal Cost Glide Trajectory")
            
            ax0[0].plot(Q_sol[0, :], Q_sol[2, :], 'k--', alpha = .5);
            ax0[0].scatter(Q_sol[0, 0], Q_sol[2, 0], marker = 'o', color = 'k', label = 'start');
            ax0[0].scatter(Q_sol[0, -1], Q_sol[2, -1], marker = 's', color = 'k', label = 'end');
            ax0[0].set_xlabel("x [m]");
            ax0[0].set_ylabel("y [m]");
            ax0[0].set_title("Centroid Trajectory");
            ax0[0].legend(loc = 'best')
            ax0[0].set_aspect('equal')
            ax0[0].axis('square')
            
            ax0[1].plot(T_VECTOR_sol, Q_sol[0, :], label = "x [m]")
            ax0[1].plot(T_VECTOR_sol, Q_sol[2, :], label = "y [m]")
            ax0[1].plot(T_VECTOR_sol, np.rad2deg(Q_sol[4, :]), label = "phi [deg]")
            ax0[1].plot(T_VECTOR_sol, Q_sol[1, :], label = "x_dot [m/s]")
            ax0[1].plot(T_VECTOR_sol, Q_sol[3, :], label = "y_dot [m/s]")
            ax0[1].plot(T_VECTOR_sol, np.rad2deg(Q_sol[5, :]), label = "phi_dot [deg/s]")
            ax0[1].set_xlabel("Time [sec]");
            ax0[1].set_ylabel("Value");
            ax0[1].set_title("State Evolution");
            ax0[1].legend(loc = 'best')
            
            
            ax0[2].plot(T_VECTOR_sol, np.rad2deg(Q_sol[6, :]), label = "theta1 [deg]")
            ax0[2].plot(T_VECTOR_sol, np.rad2deg(Q_sol[7, :]), label = "theta2 [deg]")
            ax0[2].set_xlabel("Time [sec]");
            ax0[2].set_ylabel("Value");
            ax0[2].set_title("Control Evolution");
            ax0[2].legend(loc = 'best')
        

        
    










    def animateHorizon(self, saveAsFilename = "minimal_cost_horizon_animation"):
        

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=10, bitrate = -1)
    
        # true geometry dimensions
        body_width= .30; # m
        flap_width = .15; # m
        line_thickness = 3;
        
        # scaled dimensions for drawing
        drawing_scaling = np.max([self.Q_sol[0, :], self.Q_sol[2, :]])*.15/body_width;
        body_width *= drawing_scaling;
        flap_width *= drawing_scaling;
        

        
        fig2, ax2 = plt.subplots(1, 3, figsize=(18, 6))
        
        
        
        # Helper function to get the points on the body with attitude and position changes
        def getBodyPoints(body_center_x_coord, body_center_y_coord, body_roll_angle, theta1_val, theta2_val):
            body_left_point = np.array([[-body_width/2], [0]]);
            body_right_point = np.array([[body_width/2], [0]]);
            
            flap1_tip = body_left_point + np.array([[flap_width*np.cos(np.pi-theta1_val)], [flap_width*np.sin(np.pi-theta1_val)]]);
            flap2_tip = body_right_point + np.array([[flap_width*np.cos(theta2_val)], [flap_width*np.sin(theta2_val)]]);
            
        
            rot_matrix = np.array([[np.cos(body_roll_angle), -np.sin(body_roll_angle)],
                                   [np.sin(body_roll_angle), np.cos(body_roll_angle)]])
        
            
            # note: these rotations only work because the vehicle's center is currently pinned to the origin
            body_left_point = rot_matrix @ body_left_point;
            body_right_point = rot_matrix @ body_right_point;
            flap1_tip = rot_matrix @ flap1_tip;
            flap2_tip = rot_matrix @ flap2_tip;
            
            # now, we shift the points according to our global (X, Y) position
            body_left_point += np.array([[body_center_x_coord], [body_center_y_coord]]);
            body_right_point += np.array([[body_center_x_coord], [body_center_y_coord]]);
            flap1_tip += np.array([[body_center_x_coord], [body_center_y_coord]]);
            flap2_tip += np.array([[body_center_x_coord], [body_center_y_coord]]);
            
            return body_left_point, body_right_point, flap1_tip, flap2_tip;
        
        
        
        # Initialize a configuration of patches
        
        # roll angle = 0
        body_line_left = plt.Line2D((self.Q_sol[0, 0], self.Q_sol[2, 0]), (-body_width/2, self.Q_sol[2, 0]), lw = line_thickness, color = 'black') 
        body_line_right = plt.Line2D((self.Q_sol[0, 0], self.Q_sol[2, 0]), (body_width/2, self.Q_sol[2, 0]), lw = line_thickness, color = 'black') 
        
        # flaps fully extended
        flap1_line = plt.Line2D((self.Q_sol[0, 0] - body_width/2 - flap_width, self.Q_sol[2, 0]),
                                (self.Q_sol[0, 0] - body_width/2, self.Q_sol[2, 0]), lw = line_thickness, color = 'blue');
             
        flap2_line = plt.Line2D((self.Q_sol[0, 0] + body_width/2, self.Q_sol[2, 0]),
                                (self.Q_sol[0, 0] + body_width/2 + flap_width, self.Q_sol[2, 0]), lw = line_thickness, color = 'red');                   
             
        
        state_time_line = plt.Line2D((1, 1), ((np.nanmin(self.Q_sol[0:4, :])), (np.nanmax(self.Q_sol[0:4, :]))), lw = 2,  color = 'black')
        control_time_line = plt.Line2D((1, 1), (np.rad2deg(np.nanmin(self.Q_sol[6:8, :])), np.rad2deg(np.nanmax(self.Q_sol[6:8, :]))), lw = 2, color = 'black')
        
        
        
        # Plot our state and control time histories for reference
        
        
        
        fig2.suptitle("Minimal Cost Glide Trajectory")
            
        ax2[0].plot(self.Q_sol[0, :], self.Q_sol[2, :], 'k--', alpha = .5);
        ax2[0].scatter(self.Q_sol[0, 0], self.Q_sol[2, 0], marker = 'o', color = 'k', label = 'start');
        ax2[0].scatter(self.Q_sol[0, -1], self.Q_sol[2, -1], marker = 's', color = 'k', label = 'end');
        ax2[0].set_xlabel("x [m]");
        ax2[0].set_ylabel("y [m]");
        
        xlim_val = 1.1*( np.abs(self.Q_sol[0, 0]) + body_width + flap_width)
        
        ax2[0].set_xlim(-xlim_val, xlim_val)
        ax2[0].set_title("Centroid Trajectory (figure not to scale)");
        ax2[0].legend(loc = 'upper left')
        ax2[0].set_aspect('equal')
        
        ax2[1].plot(self.T_VECTOR_sol, self.Q_sol[0, :], label = "x [m]")
        ax2[1].plot(self.T_VECTOR_sol, self.Q_sol[2, :], label = "y [m]")
        ax2[1].plot(self.T_VECTOR_sol, np.rad2deg(self.Q_sol[4, :]), label = "phi [deg]")
        ax2[1].plot(self.T_VECTOR_sol, self.Q_sol[1, :], label = "x_dot [m/s]")
        ax2[1].plot(self.T_VECTOR_sol, self.Q_sol[3, :], label = "y_dot [m/s]")
        ax2[1].plot(self.T_VECTOR_sol, np.rad2deg(self.Q_sol[5, :]), label = "phi_dot [deg/s]")
        ax2[1].set_xlabel("Time [sec]");
        ax2[1].set_ylabel("Value");
        ax2[1].set_title("State Evolution");
        ax2[1].legend(loc = 'upper right')
        
        
        ax2[2].plot(self.T_VECTOR_sol, np.rad2deg(self.Q_sol[6, :]), 'b', label = "theta1 [deg]")
        ax2[2].plot(self.T_VECTOR_sol, np.rad2deg(self.Q_sol[7, :]), 'r', label = "theta2 [deg]")
        ax2[2].set_xlabel("Time [sec]");
        ax2[2].set_ylabel("Value");
        ax2[2].set_title("Control Evolution");
        ax2[2].legend(loc = 'lower right')
        
        
        
        
        
            
            
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
            body_x_coord = self.Q_sol[0, i];
            body_y_coord = self.Q_sol[2, i];
            body_roll_angle = self.Q_sol[4, i];
            theta1_val = self.Q_sol[6, i];
            theta2_val = self.Q_sol[7, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_x_coord, body_y_coord, body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([body_left_point[0], body_x_coord])
            line.set_ydata([body_left_point[1], body_y_coord])
            return line,
        
        
        def animateBodyLineRight(i, line):
            body_x_coord = self.Q_sol[0, i];
            body_y_coord = self.Q_sol[2, i];
            body_roll_angle = self.Q_sol[4, i];
            theta1_val = self.Q_sol[6, i];
            theta2_val = self.Q_sol[7, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_x_coord, body_y_coord, body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([body_x_coord, body_right_point[0]])
            line.set_ydata([body_y_coord, body_right_point[1]])
            return line,
        
            
        def animateFlap1(i, line):
            body_x_coord = self.Q_sol[0, i];
            body_y_coord = self.Q_sol[2, i];
            body_roll_angle = self.Q_sol[4, i];
            theta1_val = self.Q_sol[6, i];
            theta2_val = self.Q_sol[7, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_x_coord, body_y_coord, body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([flap1_tip[0], body_left_point[0]])
            line.set_ydata([flap1_tip[1], body_left_point[1]])
            return line,
        
        
        def animateFlap2(i, line):
            body_x_coord = self.Q_sol[0, i];
            body_y_coord = self.Q_sol[2, i];
            body_roll_angle = self.Q_sol[4, i];
            theta1_val = self.Q_sol[6, i];
            theta2_val = self.Q_sol[7, i];
            body_left_point, body_right_point, flap1_tip, flap2_tip = getBodyPoints(body_x_coord, body_y_coord, body_roll_angle, theta1_val, theta2_val);
            
            line.set_xdata([body_right_point[0], flap2_tip[0]])
            line.set_ydata([body_right_point[1], flap2_tip[1]])
            return line,
        
        
        def animateStateTimeLine(i, line):
            line.set_xdata([i*self.DT, i*self.DT]);
        
        
        def animateControlTimeLine(i, line):
            line.set_xdata([i*self.DT, i*self.DT]);
        
        anim = animation.FuncAnimation(fig2, animationController, 
                                       init_func=init, 
                                       frames=len(self.T_VECTOR_sol), 
                                       fargs = (body_line_left, body_line_right, flap1_line, flap2_line, state_time_line, control_time_line),
                                       interval=self.DT*1000,
                                       blit=True)
        
        anim.save(saveAsFilename + '.gif', writer=writer)
        
        plt.show()    
          
        
        
        