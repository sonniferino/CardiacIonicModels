function [u,t,v,m,f,to,X]=bernus(num_steps)

t_fin = 500; %era 300
 
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
 
u =zeros(size(t));
v =zeros(size(t));
m =zeros(size(t));
f =zeros(size(t));
to=zeros(size(t));
X=zeros(size(t));
 
I_Na=zeros(size(t));
I_Ca=zeros(size(t));
I_to=zeros(size(t));
I_K=zeros(size(t));
I_K1=zeros(size(t));
I_Nab=zeros(size(t));
I_Cab=zeros(size(t));
I_NaK=zeros(size(t));
I_NaCa=zeros(size(t));



 
u(1)=-60.0;    % potential
v(1)=1; %0.998005629294136;
m(1)=0; %6.246268194901203e-04;
f(1)=1;%0.920365627223952;
to(1)=1;%0.999890987198580;
X(1)=0;%0.008526798595657;
 
 
I_Na(1)=I_Na_fun(u(1),m(1),v(1));
I_Ca(1)=I_Ca_fun(u(1),f(1));
I_to(1)=I_to_fun(u(1),to(1));
I_K(1)=I_K_fun(u(1),X(1));
I_K1(1)=I_K1_fun(u(1));
I_Nab(1)=I_Nab_fun(u(1));
I_Cab(1)=I_Cab_fun(u(1));
I_NaK(1)=I_NaK_fun(u(1));
I_NaCa(1)=I_NaCa_fun(u(1));
 

%u(2:20)=-60.0 ;
for i=2:2
   
    uo=u(i-1);
    vo=v(i-1);
    mo=m(i-1);
    fo=f(i-1);
    too=to(i-1);
    Xo=X(i-1);
    
    I_Na_o=I_Na(i-1);
    I_Ca_o=I_Ca(i-1);
    I_to_o=I_to(i-1);
    I_K_o=I_K(i-1);
    I_K1_o=I_K1(i-1);
    I_Nab_o=I_Nab(i-1);
    I_Cab_o=I_Cab(i-1);
    I_NaK_o=I_NaK(i-1);
    I_NaCa_o=I_NaCa(i-1);
    
 
    I_tot_o=  I_Na_o+ I_Ca_o+ I_to_o+I_K_o+I_K1_o+I_Nab_o+ I_Cab_o+ I_NaK_o+I_NaCa_o;
    
    L=-1/tau_v(uo);
    N=v_inf(uo)/tau_v(uo);
    %vs=exp(L*dt)*vo+1/L*(exp(L*dt)-1)*N;  %e' RL
    vs =vo + dt *v_fun(vo,uo); %e' AB1
    
    L=-alpha_m(uo)-beta_m(uo);
    N=alpha_m(uo);
   % ms=exp(L*dt)*mo+1/L*(exp(L*dt)-1)*N;  %% quella problematica ? la m
    ms =mo + dt *m_fun(mo,uo);
    
    L=-alpha_f(uo)-beta_f(uo);
    N=alpha_f(uo);
    %fs=exp(L*dt)*fo+1/L*(exp(L*dt)-1)*N;
    fs =fo + dt *f_fun(fo,uo);
    
    L=-alpha_to(uo)-beta_to(uo);
    N=alpha_to(uo);
    % tos=exp(L*dt)*too+1/L*(exp(L*dt)-1)*N;
    tos=too+ dt *to_fun(too,uo);
    
    L=-1/tau_X(uo);
    N=X_inf(uo)/tau_X(uo);
   % Xs=exp(L*dt)*Xo+1/L*(exp(L*dt)-1)*N;  %RL sulle gating
    Xs= Xo + dt *X_fun(Xo,uo);
 
    us=uo-dt*(I_tot_o); %passo esplicito
  
    
    I_Na_s=I_Na_fun( us , ms , vs );
    I_Ca_s=I_Ca_fun( us,fs);
    I_to_s=I_to_fun( us,tos );
    I_K_s=I_K_fun( us, Xs );
    I_K1_s=I_K1_fun( us );
    I_Nab_s=I_Nab_fun( us  );
    I_Cab_s=I_Cab_fun( us );
    I_NaK_s=I_NaK_fun( us );
    I_NaCa_s=I_NaCa_fun(us );
    
    I_tot_s=  I_Na_s+ I_Ca_s+ I_to_s+I_K_s+I_K1_s+I_Nab_s+ I_Cab_s+ I_NaK_s+I_NaCa_s;
    
    %%%%%%%%%%%%%%%%%   AGGIUNTO  %%%%%%%%%%%%%%%%%%%%%
    
    L=-1/tau_v( us );                     
    N=v_inf(us)/tau_v(us);
    No=v_inf(uo)/tau_v(uo);
    %vn=exp(L*dt)*vo+1/L*(exp(L*dt)-1)*N;
    vn =vo + (dt/2) *(v_fun(vo,uo)+v_fun(vs,us));  % AM2
    %vn=exp(L*dt)*vo+dt*phi1(L*dt)*No +dt^2*phi2(L*dt)*(N-No);
    
    
   % L=-alpha_m(us)-beta_m( us );
   % N=alpha_m( us );
    %No=alpha_m( uo );
    %mn=exp(L*dt)*mo+1/L*(exp(L*dt)-1)*N;  %% quella problematica ? la m
    mn =mo + (dt/2) *(m_fun(mo,uo)+m_fun(ms,us)); % AM2
    %No=alpha_m( uo );
    %mn=ms+1/L*(exp(L*dt)-1)*N;
    
    %L=-alpha_f( us )-beta_f(us );
    %N=alpha_f(us);
    %fn=exp(L*dt)*fo+1/L*(exp(L*dt)-1)*N;
    fn =fo + (dt/2) *(f_fun(fo,uo)+f_fun(fs,us));
    
    
    %L=-alpha_to( us )-beta_to( us );
    %N=alpha_to( us );
    %ton=exp(L*dt)*too+1/L*(exp(L*dt)-1)*N;
    ton =too + (dt/2) *(to_fun(too,uo)+to_fun(tos,us));
    
    
    %L=-1/tau_X( us );
    %N=X_inf( us )/tau_X( us) ;
    %Xn=exp(L*dt)*Xo+1/L*(exp(L*dt)-1)*N;
    Xn =Xo + (dt/2) *(X_fun(Xo,uo)+X_fun(Xs,us));
    
    
    
    un=uo-dt*((1/2)*I_tot_o +(1/2)*I_tot_s);
    
    u(i)=un;
    v(i)=vn;
    m(i)=mn;
    f(i)=fn;
    to(i)=ton;
    X(i)=Xn;
 
    I_Na(i)=I_Na_fun(u(i),m(i),v(i));
    I_Ca(i)=I_Ca_fun(u(i),f(i));
    I_to(i)=I_to_fun(u(i),to(i));
    I_K(i)=I_K_fun(u(i),X(i));
    I_K1(i)=I_K1_fun(u(i));
    I_Nab(i)=I_Nab_fun( u(i) );
    I_Cab(i)=I_Cab_fun( u(i) );
    I_NaK(i)=I_NaK_fun( u(i) );
    I_NaCa(i)=I_NaCa_fun( u(i) );
 
end
 
for i=3:length(t)
   
    
    uolder=u(i-2);
    volder=v(i-2);
    molder=m(i-2);
    folder=f(i-2);
    toolder=to(i-2);
    Xolder=X(i-2);
  
    uo=u(i-1);
    vo=v(i-1);
    mo=m(i-1);
    fo=f(i-1);
    too=to(i-1);
    Xo=X(i-1);
    
    I_Na_o=I_Na(i-1);
    I_Ca_o=I_Ca(i-1);
    I_to_o=I_to(i-1);
    I_K_o=I_K(i-1);
    I_K1_o=I_K1(i-1);
    I_Nab_o=I_Nab(i-1);
    I_Cab_o=I_Cab(i-1);
    I_NaK_o=I_NaK(i-1);
    I_NaCa_o=I_NaCa(i-1);
    
    I_tot_o=  I_Na_o+ I_Ca_o+ I_to_o+I_K_o+I_K1_o+I_Nab_o+ I_Cab_o+ I_NaK_o+I_NaCa_o;
    
    I_Na_older=I_Na(i-2);
    I_Ca_older=I_Ca(i-2);
    I_to_older=I_to(i-2);
    I_K_older=I_K(i-2);
    I_K1_older=I_K1(i-2);
    I_Nab_older=I_Nab(i-2);
    I_Cab_older=I_Cab(i-2);
    I_NaK_older=I_NaK(i-2);
    I_NaCa_older=I_NaCa(i-2);
    
    
    
    I_tot_older=  I_Na_older+ I_Ca_older+ I_to_older+I_K_older+I_K1_older+I_Nab_older+ I_Cab_older+ I_NaK_older+I_NaCa_older;
    
    %%%%%%%%%%%%%%%%%%%% computation of g star
    L=-1/tau_v( uo );
    N=v_inf( uo )/tau_v( uo );
    Nolder=v_inf( uolder )/tau_v( uolder );
    %vs=exp(L*dt)*vo+1/L*(exp(L*dt)-1)*N;
    vs =vo + dt *((3/2)*v_fun(vo,uo)-(1/2)*v_fun(volder,uolder));  %AB2
    %vs= exp(L*dt)*vo+((exp(L*dt)-1)/L)*N +((exp(L*dt)-1-L*dt)/(L^2))*(N-Nolder);
  
    
    L= -alpha_m( uo )-beta_m( uo );
    N=alpha_m( uo );
    Nolder=alpha_m( uolder );
    %ms=exp(L*dt)*mo+1/L*(exp(L*dt)-1)*N;  %% quella problematica ? la m
    ms =mo + dt *((3/2)*m_fun(mo,uo)-(1/2)*m_fun(molder,uolder));  
    %ms= exp(L*dt)*mo+dt*phi1(L*dt)*N +dt*phi2(L*dt)*(N-Nolder);
    
%     alpha=-alpha_m( uo )
%     bet= -beta_m( uo )
%     L=L
%     exponnnn=exp(L*dt)
% %     
%      if (-L>1e3)
% %       fprintf('colpa della m ')
% %       
% %       c=-L*dt/0.5
% %       i
% %       return
% %       
% L
%     end
% %     
    
   % L=-alpha_f(uo)-beta_f(uo);
    %N=alpha_f(uo);
    %Nolder=alpha_f(uolder);
    %fs=exp(L*dt)*fo+1/L*(exp(L*dt)-1)*N;
    %fs =fo + dt *f_fun(fo,uo);
    fs =fo + dt *((3/2)*f_fun(fo,uo)-(1/2)*f_fun(folder,uolder));
    %fs= exp(L*dt)*fo+((exp(L*dt)-1)/L)*N +((exp(L*dt)-1-L*dt)/(L^2*dt))*(N-Nolder);
  
    L=-alpha_to(uo)-beta_to(uo);
    N=alpha_to(uo);
    Nolder=alpha_to(uolder);
    %tos=exp(L*dt)*too+1/L*(exp(L*dt)-1)*N;
    %ton=too+ dt *to_fun(too,uo);
    tos =too + dt *((3/2)*to_fun(too,uo)-(1/2)*to_fun(toolder,uolder));
    %tos= exp(L*dt)*too+((exp(L*dt)-1)/L)*N +((exp(L*dt)-1-L*dt)/(L^2*dt))*(N-Nolder);
   
  
    
    %L=-1/tau_X((3/2)*uo-(1/2)*uolder);
    %N=X_inf((3/2)*uo-(1/2)*uolder)/tau_X((3/2)*uo-(1/2)*uolder);
    %Xs=exp(L*dt)*Xo+1/L*(exp(L*dt)-1)*N;
    %Xn= Xo + dt *X_fun(Xo,uo);
    Xs =Xo + dt *((3/2)*X_fun(Xo,uo)-(1/2)*X_fun(Xolder,uolder));
  
 
    %N=-(I_Na_o+I_Ca_o+I_to_o+I_K_o+I_K1_o+I_Nab_o+I_Cab_o+I_NaK_o+I_NaCa_o);
    us=uo-dt*((3/2)*I_tot_o-(1/2)*I_tot_older);
    
    %%%%%%%%%%%%%%%%%   corrector  %%%%%%%%%%%%%%%%%%%%%
    
    %L=-1/tau_v( (5/12)*us +(8/12)*uo +(-1/12)*uolder);
    %N=v_inf((5/12)*us +(8/12)*uo +(-1/12)*uolder)/tau_v((5/12)*us +(8/12)*uo +(-1/12)*uolder);
    %vn=exp(L*dt)*vo+1/L*(exp(L*dt)-1)*N;
    vn =vo + dt *((5/12)*v_fun(vs,us)+(8/12)*v_fun(vo,uo)-(1/12)*v_fun(volder,uolder));
    
    
    L=-alpha_m( uo )-beta_m( uo  );
    N=alpha_m( uo  );
    Nolder=alpha_m ( uolder );
    Nstar=alpha_m ( us );
    %mn=exp(L*dt)*mo+1/L*(exp(L*dt)-1)*N;  %% quella problematica ? la m
    %mn=exp(L*dt)*mo+dt*( phi1( L*dt )* N ) +dt^2*( phi2 (L*dt) *( (1/2) * Nstar - (1/2) * Nolder ))+dt^3*(phi3(L*dt)*(Nstar-2*N+Nolder)); 
    
   % mn=ms+dt*( phi2 (L*dt) *(  Nstar - 2*No + Nolder ));%% quella problematica ? la m
    mn =mo + dt *((5/12)*m_fun(ms,us)+(8/12)*m_fun(mo,uo)-(1/12)*m_fun(molder,uolder));
    
    %L=-alpha_f( (5/12)*us +(8/12)*uo +(-1/12)*uolder )-beta_f( (5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %N=alpha_f((5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %fn=exp(L*dt)*fo+1/L*(exp(L*dt)-1)*N;
    fn =fo + dt *((5/12)*f_fun(fs,us)+(8/12)*f_fun(fo,uo)-(1/12)*f_fun(folder,uolder));
    
    
    %L=-alpha_to( (5/12)*us +(8/12)*uo +(-1/12)*uolder )-beta_to( (5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %N=alpha_to((5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %ton=exp(L*dt)*too+1/L*(exp(L*dt)-1)*N;
    ton =too + dt *((5/12)*to_fun(tos,us)+(8/12)*to_fun(too,uo)-(1/12)*to_fun(toolder,uolder));
    
    
    %L=-1/tau_X( (5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %N=X_inf( (5/12)*us +(8/12)*uo +(-1/12)*uolder )/tau_X( (5/12)*us +(8/12)*uo +(-1/12)*uolder );
    %Xn=exp(L*dt)*Xo+1/L*(exp(L*dt)-1)*N;
    Xn =Xo + dt *((5/12)*X_fun(Xs,us)+(8/12)*X_fun(Xo,uo)-(1/12)*X_fun(Xolder,uolder)); %AM3
    
    I_Na_s=I_Na_fun( us , ms , vs );
    I_Ca_s=I_Ca_fun( us,fs);
    I_to_s=I_to_fun( us,tos );
    I_K_s=I_K_fun( us, Xs );
    I_K1_s=I_K1_fun( us );
    I_Nab_s=I_Nab_fun( us  );
    I_Cab_s=I_Cab_fun( us );
    I_NaK_s=I_NaK_fun( us );
    I_NaCa_s=I_NaCa_fun(us );
    
    I_tot_s=  I_Na_s+ I_Ca_s+ I_to_s+I_K_s+I_K1_s+I_Nab_s+ I_Cab_s+ I_NaK_s+I_NaCa_s;
    
    un=uo-dt*( (5/12)*I_tot_s +(8/12)*I_tot_o +(-1/12)*I_tot_older);
    
    u(i)=un;
    v(i)=vn;
    m(i)=mn;
    f(i)=fn;
    to(i)=ton;
    X(i)=Xn;
 
    I_Na(i)=I_Na_fun(u(i),m(i),v(i));
    I_Ca(i)=I_Ca_fun(u(i),f(i));
    I_to(i)=I_to_fun(u(i),to(i));
    I_K(i)=I_K_fun(u(i),X(i));
    I_K1(i)=I_K1_fun(u(i));
    I_Nab(i)=I_Nab_fun( u(i) );
    I_Cab(i)=I_Cab_fun( u(i) );
    I_NaK(i)=I_NaK_fun( u(i) );
    I_NaCa(i)=I_NaCa_fun( u(i) );
 
end






 