function [J_fi, J_so, J_si]=currents(v,w,s,u)

u_0=           0;
u_u=           1.55;
theta_v=       0.3;
theta_w=       0.13;
theta_v_m=    0.006;
theta_o=       0.006;
tau_v1_m=     60;
tau_v2_m=     1150;
tau_v_p=      1.4506;
tau_w1_m=      60;
tau_w2_m=      15;
k_w_m=         65;
u_w_m=         0.03;
tau_w_p=       200;
tau_fi=        0.11;
tau_o1=        400;
tau_o2=        6;
tau_so1=       30.0181;
tau_so2=       0.9957;
k_so=          2.0458;
u_so=          0.65;
tau_s1=        2.7342;
tau_s2=        16;
k_s=           2.0994;
u_s=           0.9087;
tau_si=        1.8875;
tau_w_inf=     0.07;
w_inf_star=    0.94;

tau_so =  tau_so1 + (tau_so2 - tau_so1) * (1 + tanh( k_so * (u - u_so) ) ) /2;

tau_o = (1-heaviside(u - theta_o ) )*tau_o1 +  heaviside(u - theta_o ) * tau_o2;

J_fi = -v*heaviside(u - theta_v ) * (u-theta_v) * (u_u - u) / tau_fi;

J_so = (u- u_0) * (1 - heaviside(u - theta_w ))/tau_o + heaviside(u - theta_w ) /tau_so;

J_si = -heaviside(u - theta_w )*w*s/tau_si;


end
