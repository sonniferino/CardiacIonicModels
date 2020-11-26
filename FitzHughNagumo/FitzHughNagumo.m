clear all
close all
 
t=linspace(0,1000,100001);

% parameters of simone modified

a=1.4e-3;

c2=0.1;
b=0.013/5;
d=2.30; %1

% parameters of dimension model

Vr=-85;
Vm=30;
Vd=Vm-Vr;
Vu=-57.6;

% computed parameters

mu1=a;
mu2=c2*100;
mu3=b;
mu4=d;


V(1)=Vr;
w(1)=0.0;
 
counter=1;
 
tau=t(2)-t(1);
 
f=@(t)(t<=2.0);
 
for i=2:length(t)
    
        Vold=V(i-1);
        wold=w(i-1);
        
        V_norm=(Vold-Vr)/Vd;
        
        Iapp=35*f(t(i-1));
        
        Vnew=Vold+tau*(Iapp-mu1*(Vold-Vr)*(Vold-Vm)*(Vold-Vu) -mu2*wold*(Vold-Vr));
        wnew=wold+tau*mu3*(V_norm-mu4*wold);


  
        V(i)=Vnew;
        w(i)=wnew;
        counter=counter+1;

    
end
 

subplot(1,2,1)
plot(t,V);

subplot(1,2,2)
plot(t,w);

