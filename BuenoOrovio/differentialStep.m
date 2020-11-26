function [v,w,s]=differentialStep(v,w,s,u,tau)

u_0=           0;
u_u=           1.55;
theta_v=       0.3;
theta_w=       0.13;
theta_v_m=    0.006;
theta_0=       0.006;
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

%%%%%%%%%%%%%%%%%         1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (u<theta_v_m)
    v_inf=1;
else
    v_inf=0;
end



tau_v_m=( 1-heaviside(u-theta_v_m) ) * tau_v1_m + heaviside(u-theta_v_m) * tau_v1_m;

vfun=  ( 1-heaviside(u-theta_v) ) * (v_inf - v) / tau_v_m -  heaviside(u-theta_v)*v /tau_v_p ;

%%%%%%%%%%%%%%%%%         2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w_inf = (1-heaviside(u-theta_0) )*(1 - u/tau_w_inf) + heaviside(u-theta_0) * w_inf_star;

tau_w_m=  tau_w1_m + (tau_w2_m - tau_w1_m) * (1 + tanh( k_w_m * (u - u_w_m) ) ) /2;

wfun=  ( 1-heaviside(u-theta_w) ) * (w_inf - w) / tau_w_m -  heaviside(u-theta_w)*w /tau_w_p ;

%%%%%%%%%%%%%%%%%         3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_s = (1-heaviside(u - theta_w ) )*tau_s1 +  heaviside(u - theta_w ) * tau_s2;

sfun= (( 1 + tanh(k_s *( u-u_s)) )/2-s)/tau_s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=v+tau*vfun;
w=w+tau*wfun;
s=s+tau*sfun;

end