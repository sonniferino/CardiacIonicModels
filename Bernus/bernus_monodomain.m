function [u,t]=bernus_monodomain(num_steps,t_fin)

%%% VERSIONE CORRETTA
%%% LA VELOCITA DI CONDUZIONE ? circa 0.77 m/s
%%% NOTA la conduttivita deve essere in mS/cm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   parameters
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
T=273.15+37;
R=8314.4598;  % J/(mol*K)
F=96485.33289; % C/mol
 
Ca_i=0.0004;
Ca_e=2.0;
Na_i=10;
Na_e=138;
K_i=140;
K_e=4;
 
g_Na=16;
g_Ca=0.064;
g_to=0.4;
g_K=0.019;
g_K1=3.9;
g_Nab=0.001;
g_Cab=0.00085;
g_NaK=1.3;
g_NaCa=1000;
 
E_Na = R*T/F*log(Na_e/Na_i);
E_Ca = R*T/2/F*log(Ca_e/Ca_i);
E_to = R*T/F*log( ( 0.043*Na_e+K_e)/(0.043*Na_i+K_i) );
E_K  = R*T/F*log(K_e/K_i);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Na_fun
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%   functions for v
 
v_inf=@(u)( 0.5*(1.0-tanh(7.74+0.12*u))  );
tau_v=@(u)( 0.25 +2.24*(( 1-tanh(7.74+0.12*u)) ./ ( 1-tanh(0.07*(u+92.4)) )) );
v_fun=@(v,u)( (v_inf(u)-v)./tau_v(u) );
 
%%%%%%%%%%   functions for m
 
alpha_m=@(u)( (0.32*(u+47.13))./( 1-exp(-0.1*(u+47.13) ) ) );
beta_m=@(u)( 0.08*exp(-u/11) );
m_fun=@(m,u)(alpha_m(u)*(1-m)- beta_m(u)*m );
 
%%%%%%%%%%   I_Na_fun
 
I_Na_fun=@(u,m,v)(g_Na*m^3*v^2*(u-E_Na));
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Ca_fun
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
f_Ca=1/(1+Ca_i/0.0006);
 
 
alpha_d=@(u)( (14.98*exp(-0.5*((u-22.36)/16.68)^2))/(16.68*sqrt(2*pi)) );
beta_d=@(u)( 0.1471-(5.3*exp(-0.5*((u-6.27)/14.93)^2))/(14.93*sqrt(2*pi)) );
 
d_inf=@(u)( alpha_d(u)/( alpha_d(u) + beta_d(u)  ) );
 
alpha_f=@(u)( (6.87*1e-3)/(1+exp(-(6.1546-u)/6.12)) );
beta_f=@(u)( (0.069*exp(-0.11*(u+9.825))+0.011)/(1+exp(-0.278*(u+9.825)))+5.75*1e-4 );
 
f_fun=@(f,u)(alpha_f(u)*(1-f)- beta_f(u)*f );
 
I_Ca_fun=@(u,f)( g_Ca*f_Ca * d_inf(u)*f *(u-E_Ca) );
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   To
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
alpha_r=@(u)( (0.5266*exp(-0.0167*(u-42.2912)))/(1+exp(-0.0943*(u-42.2912))) );
beta_r=@(u)( (5.186*1e-5*u+0.5149*exp(-0.1344*(u-5.0027)))/(1+exp(-0.1348*(u-5.186*1e-5))) );
 
r_inf=@(u)(alpha_r(u)/(alpha_r(u) + beta_r(u) ) );
 
alpha_to=@(u)( (5.612e-5*u+0.0721*exp(-0.173*(u+34.2531)))/(1+exp(-0.1732*(u+34.2531))) );
beta_to=@(u)(  (1.215e-4*u+0.0767*exp(-1.66e-9*(u+34.0235)))/(1+exp(-0.1604*(u+34.0235))) );
 
to_fun=@(to,u)( alpha_to(u) * (1-to) - beta_to(u) * to);
 
I_to_fun=@(u,to)(g_to*r_inf(u)*to*(u-E_to));
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_K
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
tau_prime_X=@(u)(40*(1-tanh( 160+2*u ) ));
tau_X=@(u)( 240*exp(-(25.5+u)^2/156)+182*(1+tanh(0.154+0.0116*u))+tau_prime_X(u) );
X_inf=@(u) ( 0.988/(1+exp(-0.861-0.0620*u)) );
 
X_fun=@(X,u)( ( X_inf(u)-X )/tau_X(u) );
 
I_K_fun=@(u,X)( g_K*X^2*(u-E_K) );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_K1
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
alpha_K1=@(u)(0.1/(1+exp(0.06*(u-E_K-200))) );
beta_K1=@(u)( (3*exp(2*1e-4*(u-E_K+100)) + exp(0.1*(u-E_K-10)))/(1+exp(-0.5*(u-E_K))) );
K1_inf=@(u) ( alpha_K1(u)/(alpha_K1(u) + beta_K1(u) ));
 
I_K1_fun=@(u)( g_K1*K1_inf(u)*(u-E_K) );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Na,b
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
I_Nab_fun=@(u)(g_Nab*(u-E_Na));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Ca,b
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
I_Cab_fun=@(u)(g_Cab*(u-E_Ca));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Na,K
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
sigma=0.1428*(exp(Na_e/67.3)-1);
 
f_NaK=@(u)( 1/(1+0.1245*exp(-0.0037*u)+0.0365*sigma*exp(-0.037*u)) );
 
f_prime_NaK=@(u)( (1/(1+(10/Na_i)^(1.5) ))*( K_e/(K_e +1.5)) );
 
I_NaK_fun=@(u)( g_NaK*f_NaK(u)*f_prime_NaK(u) );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   I_Na,Ca
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
f_NaCa=@(u)( (87.5^3+Na_e^3)^(-1)*(1.38+Ca_e)^(-1)*(1+0.1*exp(-0.024*u))^(-1)*(Na_i^3*Ca_e*exp(0.013*u)-Na_e^3*Ca_i*exp(-0.024*u)) );
I_NaCa_fun=@(u) (g_NaCa * f_NaCa(u) );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%   init simulation
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
t=linspace(0,t_fin,t_fin*num_steps+1);
dt=t(2)-t(1);
x=linspace(0,1,400+1);
u=zeros(size(x))';
u(1:end)=-86.2;
sss=length(u)-1;
u(1:sss*0.04+1)=30;

  

h=diff(x)';
M=1/3*diag([h;0]+[0;h])+1/6*diag(h,1)+1/6*diag(h,-1);
A=diag([1./h;0]+[0;1./h])-diag(1./h,1)-diag(1./h,-1);

Cm=2.0;             %muF / cm^2
chi=2000;           %cm^-1
sigma=12*10.3429/(12.0+10.3429); %5.6; % 


% chi=50;
% sigma=0.55;
% Cm=0.02;  
% sigma/(chi*Cm)
%  break
  


% 
%  for i=1:length(t)
%      
%     quale=find(tref==t(i));
%     u(i)=Vref(quale);
%  end
% 

v =zeros(size(x))';
m =zeros(size(x))';
f =zeros(size(x))';
to=zeros(size(x))';
X=zeros(size(x))';
 
I_Na=zeros(size(x))';
I_Ca=zeros(size(x))';
I_to=zeros(size(x))';
I_K=zeros(size(x))';
I_K1=zeros(size(x))';
I_Nab=zeros(size(x))';
I_Cab=zeros(size(x))';
I_NaK=zeros(size(x))';
I_NaCa=zeros(size(x))';


Iion=zeros(size(x))';
 
v(:)=1;
m(:)=0;
f(:)=1;
to(:)=1;
X(:)=0;


for i=2:length(t)

    for j=1:length(x)
   
    uo=u(j);
    vo=v(j);
    mo=m(j);
    fo=f(j);
    too=to(j);
    Xo=X(j);
    
    I_Na_o=I_Na_fun(uo,mo,vo);
    I_Ca_o=I_Ca_fun(uo,fo);
    I_to_o=I_to_fun(uo,too);
    I_K_o=I_K_fun(uo,Xo);
    I_K1_o=I_K1_fun(uo);
    I_Nab_o=I_Nab_fun( uo );
    I_Cab_o=I_Cab_fun( uo );
    I_NaK_o=I_NaK_fun( uo );
    I_NaCa_o=I_NaCa_fun( uo );
    
    L=-1/tau_v(uo);
    N=v_inf(uo)/tau_v(uo);
    vn=exp(L*dt)*vo+1/L*(exp(L*dt)-1)*N;
    
    L=-alpha_m(uo)-beta_m(uo);
    N=alpha_m(uo);
    mn=exp(L*dt)*mo+1/L*(exp(L*dt)-1)*N;  %% quella problematica ? la m
    
    
    L=-alpha_f(uo)-beta_f(uo);
    N=alpha_f(uo);
    fn=exp(L*dt)*fo+1/L*(exp(L*dt)-1)*N;
    
    L=-alpha_to(uo)-beta_to(uo);
    N=alpha_to(uo);
    ton=exp(L*dt)*too+1/L*(exp(L*dt)-1)*N;
    
    L=-1/tau_X(uo);
    N=X_inf(uo)/tau_X(uo);
    Xn=exp(L*dt)*Xo+1/L*(exp(L*dt)-1)*N;
    
    
    v(j)=vn;
    m(j)=mn;
    f(j)=fn;
    to(j)=ton;
    X(j)=Xn;
 
    
    
   Iion(j)=-1.0*(I_Na_o+I_Ca_o+I_to_o+I_K_o+I_K1_o+I_Nab_o+I_Cab_o+I_NaK_o+I_NaCa_o);
 
    end

     %u=  (M+dt*0.01*A)\(M*u+dt*M*Iion);
    % u=(M+dt*0.01*A)\(M*u+dt*M*Iion);
     
    u =(chi*Cm*M+dt*sigma*A)\(chi*Cm*M*u +dt*chi*Cm*M*Iion); % IMPLICIT EUL
    %u =(chi*Cm*M+0.5*dt*sigma*A)\(chi*Cm*M*u -0.5*dt*sigma*A*u +dt*chi*Cm*M*Iion);  %% CRANK NICOLSON
     
    plot(x,u)
    ylim([-90 60])
    title(num2str((i-1)*dt))
    pause(0.000001)
    
    
end


end