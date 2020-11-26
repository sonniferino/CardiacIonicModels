function [Ta,t]=active_stress_function(num_steps,u)
t_fin = 500; %era 300
t=linspace(0,t_fin,t_fin*num_steps+1);
dt=t(2)-t(1);


V =zeros(size(t));
Ta=zeros(size(t));
Ta(1)=0;
%%%%%%%%%%%% active stress computation
Vrest=-90.272;
kTa=0.6759;
Vmax=42;
%%% con nash funzionano
 eps0=1;
 epsinf=0.01;


nash=0;
rho=0.3;




for i=2:length(t)
    
V(i)=(u(i) - Vrest) / (Vmax - Vrest);
    
if(nash==1)
    if(V(i)<0.05)
    eps=eps0;
    else
        eps=epsinf;
    end
else
    eps=eps0+(epsinf-eps0)*exp(-exp(-rho*(u(i)-Vrest)));
end


Ta(i) = (dt * eps * kTa * (u(i)-Vrest) + Ta(i-1)) / (1.0 + dt * eps);

%Ta(i) = (dt * eps * kTa * V(i) + Ta(i-1)) / (1.0 + dt * eps);


end
%close all;

max(Ta)
 plot(t,u, 'r');
hold on 
plot(t,Ta, 'r');
