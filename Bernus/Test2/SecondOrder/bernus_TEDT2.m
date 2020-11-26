function [u,t,v,m,f,to,X]=bernus_TEDT2(num_steps,tref,Vref)  
%,I_Na,I_Ca,I_to, I_K,I_K1,I_Nab,I_Cab,I_NaK,I_NaCa
%,uref,tref, I_Naref,I_Caref,I_toref, I_Kref,I_K1ref,I_Nabref,I_Cabref,I_NaKref,I_NaCaref
%

t_fin = 400;
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
tau=t(2)-t(1);
 
u =zeros(size(t));
v =zeros(size(t));
m =zeros(size(t));
f =zeros(size(t));
to=zeros(size(t));
X=zeros(size(t));


 
u(1)=-60.0;    % potential
v(1)=1; %0.998005629294136;
m(1)=0; %6.246268194901203e-04;
f(1)=1;%0.920365627223952;
to(1)=1;%0.999890987198580;
X(1)=0;%0.008526798595657;


ratio=(length(tref)-1)/(length(t)-1);   
u=Vref(1:ratio:end);
 
for i=2:length(t)
    
    un=u(i);
    uo=u(i-1);
    vo=v(i-1);
    mo=m(i-1);
    fo=f(i-1);
    too=to(i-1);
    Xo=X(i-1);

    
    %indo=(i-2)*ratio+1;
    %indn=(i-1)*ratio+1;
    %indm=(2*i-3)*ratio/2 +1;
    
   
    L=-1/tau_v( un );
    Lo=-1/tau_v( uo );

    N=v_inf(un)/tau_v(un);
    No=v_inf(uo)/tau_v(uo);
   
    alpha=(L-Lo)/tau;
    beta=Lo-t(i-1)*(L-Lo)/tau;
    integraleLp=@(t)(alpha/2)*t^2+beta*t;
  
    Lp= integraleLp(t(i))-integraleLp(t(i-1));  
     
    fu=@(w) exp(Lo*(t(i)-w)-((L-Lo)/(2*tau))*(w-t(i-1)).^2);
    fu_der=@(w) (-L-((L-Lo)/(tau)).*(w-t(i-1))) * fu ( w ) ;
    
    g=@(w) No+(N-No)*(w-t(i-1))/tau;
    g_der=@(t)(N-No)/tau;
    
    cnst=exp(tau*(L-Lo)./2);
    
    funzione=@(w)cnst*fu(w)*g(w);
    der_funzione=@(w)cnst*(fu(w)*g_der(w)+fu_der(w)*g(w));
     
    integrale=(tau/2)*(funzione(t(i))+funzione(t(i-1)))+(tau^2/12)*(der_funzione(t(i-1))- der_funzione(t(i)));
     
    % esplicito
    vn=exp(Lp)*vo+integrale;
    %%%%%%%%%%%%%%%%%%



    L=-alpha_m(un)-beta_m( un );
    Lo=-alpha_m(uo)-beta_m( uo );
   
    N=alpha_m( un );
    No=alpha_m( uo );
    
    alpha=(L-Lo)/tau;
    beta=Lo-t(i-1)*(L-Lo)/tau;
    integraleLp=@(t)(alpha/2)*t^2+beta*t;
  
    Lp= integraleLp(t(i))-integraleLp(t(i-1));  
     
    fu=@(w) exp(Lo*(t(i)-w)-((L-Lo)/(2*tau))*(w-t(i-1)).^2);
    fu_der=@(w) (-L-((L-Lo)/(tau)).*(w-t(i-1))) * fu ( w ) ;
    
    g=@(w) No+(N-No)*(w-t(i-1))/tau;
    g_der=@(t)(N-No)/tau;
    
    cnst=exp(tau*(L-Lo)./2);
    
    funzione=@(w)cnst*fu(w)*g(w);
    der_funzione=@(w)cnst*(fu(w)*g_der(w)+fu_der(w)*g(w));
     
    integrale=(tau/2)*(funzione(t(i))+funzione(t(i-1)))+(tau^2/12)*(der_funzione(t(i-1))- der_funzione(t(i)));
     
    % esplicito
    mn=exp(Lp)*mo+integrale;
   
    L=-alpha_f( un )-beta_f(un );
    Lo=-alpha_f( uo )-beta_f(uo );
    
    N=alpha_f(un);
    No=alpha_f(uo);
 
    alpha=(L-Lo)/tau;
    beta=Lo-t(i-1)*(L-Lo)/tau;
    integraleLp=@(t)(alpha/2)*t^2+beta*t;
  
    Lp= integraleLp(t(i))-integraleLp(t(i-1));  
     
    fu=@(w) exp(Lo*(t(i)-w)-((L-Lo)/(2*tau))*(w-t(i-1)).^2);
    fu_der=@(w) (-L-((L-Lo)/(tau)).*(w-t(i-1))) * fu ( w ) ;
    
    g=@(w) No+(N-No)*(w-t(i-1))/tau;
    g_der=@(t)(N-No)/tau;
    
    cnst=exp(tau*(L-Lo)./2);
    
    funzione=@(w)cnst*fu(w)*g(w);
    der_funzione=@(w)cnst*(fu(w)*g_der(w)+fu_der(w)*g(w));
     
    integrale=(tau/2)*(funzione(t(i))+funzione(t(i-1)))+(tau^2/12)*(der_funzione(t(i-1))- der_funzione(t(i)));
     
    % esplicito
    fn=exp(Lp)*fo+integrale;
    
    
    L=-alpha_to( un )-beta_to( un );
    Lo=-alpha_to( uo )-beta_to( uo );
    
    N=alpha_to( un );
    No=alpha_to( uo );
    
    alpha=(L-Lo)/tau;
    beta=Lo-t(i-1)*(L-Lo)/tau;
    integraleLp=@(t)(alpha/2)*t^2+beta*t;
  
    Lp= integraleLp(t(i))-integraleLp(t(i-1));  
     
    fu=@(w) exp(Lo*(t(i)-w)-((L-Lo)/(2*tau))*(w-t(i-1)).^2);
    fu_der=@(w) (-L-((L-Lo)/(tau)).*(w-t(i-1))) * fu ( w ) ;
    
    g=@(w) No+(N-No)*(w-t(i-1))/tau;
    g_der=@(t)(N-No)/tau;
    
    cnst=exp(tau*(L-Lo)./2);
    
    funzione=@(w)cnst*fu(w)*g(w);
    der_funzione=@(w)cnst*(fu(w)*g_der(w)+fu_der(w)*g(w));
     
    integrale=(tau/2)*(funzione(t(i))+funzione(t(i-1)))+(tau^2/12)*(der_funzione(t(i-1))- der_funzione(t(i)));
     
    % esplicito
    ton=exp(Lp)*too+integrale;
    
    L=-1/tau_X( un );
    Lo=-1/tau_X( uo );
    
    N=X_inf( un )/tau_X( un) ;
    No=X_inf( uo )/tau_X( uo) ;
    
   
    alpha=(L-Lo)/tau;
    beta=Lo-t(i-1)*(L-Lo)/tau;
    integraleLp=@(t)(alpha/2)*t^2+beta*t;
  
    Lp= integraleLp(t(i))-integraleLp(t(i-1));  
     
    fu=@(w) exp(Lo*(t(i)-w)-((L-Lo)/(2*tau))*(w-t(i-1)).^2);
    fu_der=@(w) (-L-((L-Lo)/(tau)).*(w-t(i-1))) * fu ( w ) ;
    
    g=@(w) No+(N-No)*(w-t(i-1))/tau;
    g_der=@(t)(N-No)/tau;
    
    cnst=exp(tau*(L-Lo)./2);
    
    funzione=@(w)cnst*fu(w)*g(w);
    der_funzione=@(w)cnst*(fu(w)*g_der(w)+fu_der(w)*g(w));
     
    integrale=(tau/2)*(funzione(t(i))+funzione(t(i-1)))+(tau^2/12)*(der_funzione(t(i-1))- der_funzione(t(i)));
     
    % esplicito
    Xn=exp(Lp)*Xo+integrale;
 
    v(i)=vn;
    m(i)=mn;
    f(i)=fn;
    to(i)=ton;
    X(i)=Xn;
  
end
 