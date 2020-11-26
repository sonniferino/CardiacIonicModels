clear all
close all
 
t=linspace(0,1000,100001);

% parameters of dimensionless model

c1=0.26;
c2=0.1;
a=0.13;
b=0.013;
d=1;

% parameters of dimension model

Vr=-80;
Vm=30;
Vd=Vm-Vr;
Vu=(1-a)*Vr+a*Vm;

% computed parameters

mu1=c1/Vd^2;
mu2=c2;
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
        
        Iapp=30*f(t(i-1));
        
        Vnew=Vold+tau*(Iapp-mu1*(Vold-Vr)*(Vold-Vm)*(Vold-Vu) -mu2*wold*(Vold-Vr));
        wnew=wold+tau*mu3*(V_norm-mu4*wold);


  
        V(i)=Vnew;
        w(i)=wnew;
        counter=counter+1;

    
end
 

hold on 
plot(t,V);

figure(2)
plot(t,w);

