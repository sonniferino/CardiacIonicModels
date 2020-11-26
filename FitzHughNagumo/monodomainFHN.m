clear all
close all




%%%%%%%%% param pezzuto %
sigmai=1.7;
sigmae=6.2;

sigma=(sigmai*sigmae)/(sigmai+sigmae);
Cm = 1.0;
Vr = -85;
Vm = 30;
Vu = -57.6;
chi=1400;
a=1.4e-3;

Vd=Vm-Vr;


%%%%%%%%%%
t=linspace(0,1000,10001);
tau=t(2)-t(1);

x=linspace(0,2,2001);
u=zeros(size(x))';
u(1:end)=Vr;
%u(1:50)=30;

uo=ones(size(x))';
Iion=zeros(size(x))';
Iapp=zeros(size(x))';


%f=@(x,t)200*(t<=2)*(x<=0.2);
f=@(x,t)(t<=1.98)*(x<=0.2)*180.0*exp(1/((t-1.0)*(t-1.0)-1.0));

h=diff(x)';
M=1/3*diag([h;0]+[0;h])+1/6*diag(h,1)+1/6*diag(h,-1);
A=diag([1./h;0]+[0;1./h])-diag(1./h,1)-diag(1./h,-1);

A=(tau*sigma/chi)*A;

for i=2:length(t)
    
    t(i)
    
    for j=1:length(x)
                 
        uo(j)=u(j);
        
        if (t(i)<=1.98)
        Iapp(j)=f(x(j),t(i));
        else
            Iapp(j)=0;
        end
          

       Iion(j)=a*(uo(j)-Vr)*(uo(j)-Vu)*(uo(j)-Vm);
    
    end
    
    

    MAT=Cm*M+A;
    rhs=-tau*M*Iion+Cm*M*u+(tau/chi)*Iapp;
    
    u=MAT\rhs;
   
    plot(x,u) 
    pause(0.01)
    
%     if t(i)>9
%       pause(1)
%     end
end

