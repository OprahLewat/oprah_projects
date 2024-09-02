# Simple model with pure counterflow 
# Written by Oprah Lewat
# 11 August 2024
import numpy as np 
import math 
import scipy.optimize
import matplotlib.pyplot as plt 
import CoolProp 
def solve_Qdot_dp (N, L, output=0, debug=0):
    #Geometry
    d_i = 29.1e-3
    d_o = 33.4e-3
    A_i = math.pi*d_i*L
    A_o = math.pi*d_o*L
    A = (math.pi/4)*(d_i**2)
    L_tot = L*(N**2)
    #Waterside
    rho_w= 998.3 #kg/me3
    mu_w = 1.003e-3 #kg/ms
    cp_w = 4184 #J/kgK
    k_w = 0.598 #W/mK
    T_wi = 20 #C
    v_w = 1 #m/s
    eps = 58e-6 #m(epsilon)
    m_dot = rho_w*A*v_w
    Re_w = (rho_w*v_w*d_i)/mu_w
    Pr_w = (mu_w*cp_w)/k_w
    eps_r= eps/d_i
    f = 0.25*(math.log10(0.27*eps_r+(5.74/(Re_w)**.9))**-2)
    Nu_w = ((f/8)*(Re_w-1000)*Pr_w)/(1+12.7*((f/8)**0.5)*(((Pr_w)**(2/3))-1))
    h_i = (Nu_w*k_w)/d_i
    #Tube Wall
    k_wall = 15 #W/mK
    r_ratio = d_o/d_i
    R_wall = (math.log(r_ratio))/(2*math.pi*k_wall*L) #K/w
    #Steam side
    T_o = 40.29 #C
    vf = 0.001008 #(m*e3)/kg
    vg = 19.233 #m*e3/kg
    rho_f = 1/vf #kg/m*e3
    rho_g = 1/vg #kg/m*e3
    h_fg = 2485300 #J/kg

    #Binary search to find outside wall temperature T_so
    T_so_min = T_wi
    T_so_max = T_o
    err = 1e10
    it = 0
    while abs(err) > 1e-4 and it<=20 :
        it = it+1
        T_so = 0.5*(T_so_min+T_so_max)
        # Outside convective heat transfer coefficient 
        Ja = (cp_w*(T_o-T_so))/h_fg
        h_fg_s = h_fg*(1+(0.68*Ja))
        delta_rho = rho_f - rho_g
        delta_T = T_o - T_so
        g = 9.81
        h_o = 0.729*(((g*rho_f*delta_rho*h_fg_s*((k_w**3))/mu_w*delta_T*d_o*N))**0.25)
       
       #Overall heat transfer
        UA = ((1/(h_o*A_o))+R_wall+(1+(h_i*A_i)))**-1
        C_min = rho_w*v_w*A*cp_w
        NTU = UA/ C_min
        Q_max = C_min*(T_o-T_wi)
        epsilon = (1-(math.exp(-NTU)))
        Q_dot = epsilon*Q_max*N*N

        #Heat transfer through water film
        Q_dot_film = h_o*A_o*delta_T*N*N

        #Water outlet temperature 
        UA_w = (R_wall+(1/(h_i*A_i)))**(-1)
        Q_wdot = UA_w*(T_so-T_wi)
        T_we = (Q_wdot/(m_dot*cp_w))+T_wi

        #Find normalized error 
        err = (Q_dot-Q_dot_film)/Q_dot
        if err > 0 :
            T_so_max = T_so
        else:
            T_so_min = T_so

        #Pressure drop
        dp_0 = f*(L/d_i)*0.5*rho_w*(v_w**2)

        #Print iteration information 
        if debug == 1:
            print('{:3.0f}'.format(it),'{:10.3e}'.format(err),'{:10.3f}[\u00B0C]'.format(T_so),'|{:10.3f}[\u00B0C]'.format(T_so_max))

        if output == 1:
            print()
            print ('N = {:.0f}'.format(N), 'L = {:.3f}[m]'.format(L))
            print('Re = {:.0f}'.format(Re_w), 'Pr = {:.0f}'.format((Pr_w)))
            print('UA = {:.3f}[W/mK]'.format(UA), 'h_i = {:.3f}[W/m\u00b2K]'.format(h_i), 'h_o = {:.3f}[W/m\u00b2K]'.format(h_o))
            print('T_wi = {:.1f} [\u00B0C]'.format(T_wi), 'T_we = {:.1f} [\u00B0C]'.format(T_we))
            print('T_sat = {:.1f} [\u00B0C]'.format(T_o), 'T_so = {:.1f} [\u00B0C]'.format(T_so))
            print("epsilon = {:.4f}".format(epsilon), 'Q_dot_film = {:.3f}[kW]'.format(Q_dot_film/1e3), 'Q_w_dot = {:.3f}[kW]'.format(Q_wdot))
            print('Q_dot = {:.3f} [kW]'.format(Q_dot/1e3), 'dp_0 = {:.3f} [kPa]'.format(dp_0/1e3), 'L_tot = {:.3f}[m]'.format(L_tot))

        return Q_dot, dp_0, L_tot, L
solve_Qdot_dp(4,1.692,1,1)
