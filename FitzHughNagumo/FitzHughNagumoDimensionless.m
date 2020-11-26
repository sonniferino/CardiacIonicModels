


clear all
close all
 
t=linspace(0,1000,100000);
 
V(1)=0;
w(1)=0.0;
 
counter=1;
 
tau=t(2)-t(1);
Iapp=0;


a = 0.13;
mu1 = 0.26;
mu2 = 0.1 ;
mu3 = 0.013;
mu4 = 1;

f=@(t)( t<=2.0 );
 
for i=2:length(t)
 
    Vold=V(i-1);
    wold=w(i-1);
    
    Iapp=0.350*f(t(i-1));
      
    wnew= wold+ tau*mu3*(Vold-mu4*wold);
    
    Vnew=Vold+tau*(Iapp-mu1*Vold*(Vold-a)*(Vold-1)-mu2*Vold* wnew);
    

    
    V(i)=Vnew;
    w(i)=wnew;
    
end
 
plot(t,V);

figure(2)
plot(t,w);