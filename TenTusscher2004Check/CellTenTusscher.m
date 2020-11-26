function [t,V_vect,Ca_i_vect,Ca_sr_vect,Na_i_vect,K_i_vect]=CellTenTusscher(num_steps,t_fin)


epicardium=1;
endo=0;
Mcell=0;

initialize_variables_and_parameters

t=linspace(0,t_fin,t_fin*num_steps+1);
tau=t(2)-t(1);

%L_vect=zeros(size(t))';

% Potential
V_vect=zeros(size(t))';
%Gating Variables
m_vect=zeros(size(t))';
h_vect=zeros(size(t))';
j_vect=zeros(size(t))';
d_vect=zeros(size(t))';
f_vect=zeros(size(t))';
fca_vect=zeros(size(t))';
r_vect=zeros(size(t))';
s_vect=zeros(size(t))';
chis_vect=zeros(size(t))';
chir1_vect=zeros(size(t))';
chir2_vect=zeros(size(t))';
g_vect=zeros(size(t))';
%Ionic Concentrations
Ca_i_vect    =zeros(size(t))';
Ca_sr_vect   =zeros(size(t))';
Na_i_vect    =zeros(size(t))';
K_i_vect     =zeros(size(t))';
Ca_i_tot_vect=zeros(size(t))';
Ca_sr_tot_vect=zeros(size(t))';


% Initializations of variables from code TT2004
V_vect(1)= -12; %-86.2; % -86.4938;
%V_vect(1)=-12;%V_vect(1)=20.0;

m_vect(1)=0;    %controllato: loro codice
h_vect(1)=0.75;   %controllato: loro codice
j_vect(1)=0.75;   %controllato: loro codice
d_vect(1)=0.0;   %controllato: loro codice
f_vect(1)=1.0;  %controllato: loro codice
fca_vect(1)=1.0;    %controllato: loro codice %attenzione va sopra l'uno
r_vect(1)=0.0;   %controllato: loro codice
s_vect(1)=1; %controllato: loro codice
chis_vect(1)=0; %controllato: loro codice
chir1_vect(1)=0;   %controllato: loro codice
chir2_vect(1)=1;  %controllato:  diverso
g_vect(1)=1.0; %controllato: loro codice

Ca_i_vect(1) =0.000080; %controllato: loro codice
Ca_sr_vect(1)=0.56;   %controllato: loro codice
Na_i_vect(1) =11.6;  %controllato: loro codice
K_i_vect(1)  =138.3;

Ca_i_bufc  =Ca_i_vect(1)*Buf_c  /(Ca_i_vect(1)+K_bufc);
Ca_sr_bufsr=Ca_sr_vect(1)*Buf_sr/(Ca_sr_vect(1)+K_bufsr);

Ca_i_tot_vect(1)=Ca_i_vect(1)+Ca_i_bufc;
Ca_sr_tot_vect(1)=Ca_sr_vect(1)+Ca_sr_bufsr;

%%% reversal potential

%tenTusscher paper
if 0
    
    
    % I valori sono un po' errati ma le unit? tornano con wiki
    R=8.3143; %J/(K*mol)
    T=310.0;  %K
    F=96.4867; % C/mmol
    
    %in questo caso RT/F esce in milliVolt
    
    % R*T/F = J K mmol / ( K mol C) = mV 
    
end
%tenTusscher code
if 1
    % I valori sono un po' errati ma le unit? tornano con wiki
    R=8314.472;     % mJ/(K*mol)
    T=310.0;        % K
    F=96485.3415;   % C/mol
    %in questo caso RT/F esce in milliVolt
end

%bernus & wiki
if (0)
    R=8.3144598;  % J/(K*mol)
    T=273.15+37;
    F=96.48533289; % C/mmol
    %in questo caso RT/F esce in  milliVolt
end

threshold = -37; % code
%threshold = -60; % paper


rhoKNa=0.03;

E_Na=@(Na) (R*T/F)*log(Na_o/Na);

E_K= @(K) (R*T/F)*log(K_o/K);

E_Ca= @(Ca)(R*T/(2*F))*log(Ca_o/Ca);

E_Ks= @(Na,K) (R*T/F)*log((K_o+rhoKNa*Na_o)/(K+rhoKNa*Na));

%%%%%%%%%%%%%%%%%%%%%% chiK1 %%%%%%%%%%%%%%%%%%%   %la prima non  integrata

alpha_chiK1=@(V,K) (    0.1/(1+exp(0.06*(V-E_K(K)-200.0)))    );
beta_chiK1=@(V,K) (    (3.0*exp(0.0002*(V-E_K(K)+100))+exp(0.1*(V-E_K(K)-10)) )/(1+exp(-0.5*(V-E_K(K))))   );

chiK1_inf=@(V,K) (   alpha_chiK1(V,K)/(alpha_chiK1(V,K)+beta_chiK1(V,K))   );


%%%%%%%%%%%%%%%%% currents %%%

I_Na= @(V,m,h,j,Na) (   G_Na*m^3*h*j*(V - E_Na(Na))   );

I_CaL=@(V,d,f,fca,Ca_i) (   G_CaL*d*f*fca*4* ((V*(F^2))/(R*T)) * ((Ca_i*exp((2*V*F)/(R*T))-0.341*Ca_o)/(exp((2*V*F)/(R*T))-1.0))   );

I_to=@(V,r,s,K) (   G_to *r*s*(V-E_K(K))   );

I_Ks=@(V,chis,Na,K)(   G_Ks * (chis)^2*(V-E_Ks(Na,K))   );

I_Kr=@(V,chir1,chir2,K)(   G_Kr* sqrt(K_o/5.4) *chir1*chir2*(V-E_K(K))   );



I_K1=@(V,K)(   G_K1 * sqrt(K_o/5.4) * chiK1_inf(V,K) *(V-E_K(K))   );


I_NaCa=@(V,Na,Ca)(   k_NaCa *(exp((gamma*V*F)/(R*T))*Na^3*Ca_o - exp((gamma-1)*(V*F)/(R*T))*Na_o^3*Ca*alpha)/((K_mNai^3+Na_o^3)*(K_mCa+Ca_o)*(1+k_sat*exp((gamma-1)*((V*F)/(R*T)))))   );

I_NaK=@(V,Na)(   P_NaK*(K_o*Na)/((K_o+K_mK)*(Na+K_mNa)*(1+0.1245*exp(-0.1*(V*F)/(R*T))+0.0353*exp(-V*F/(R*T))))    );

I_pCa=@(Ca) (   G_pCa *Ca/(K_pCa+Ca)   );

I_pK=@(V,K)(   G_pK*(V-E_K(K))/(1+exp((25.0-V)/5.98))   );

I_bNa=@(V,Na)(   G_bNa * (V-E_Na(Na))   );

I_bCa=@(V,Ca) (   G_bCa * (V-E_Ca(Ca))   );



I_leak=@(Ca,Ca_sr)(   V_leak*(Ca_sr-Ca)   );

I_up=@(Ca)(   V_maxup/(1+(K_up^2)/(Ca^2))   );

I_rel=@(d,g,Ca_sr)(  (a_rel*(Ca_sr^2)/(b_rel^2+Ca_sr^2)+c_rel)*d*g   );


%%%%%%%%%%%%%% gating %%%%%%%%%%%%

%%%%%%%%%%%%%% m %%%%%%%%%%%%%%%%

m_inf=@(V)( 1/((1+exp((-56.86-V)/9.03))^2)   );
alpha_m=@(V)( 1/(1+exp((-60-V)/5))   );
beta_m=@(V)(0.1/(1+exp((V+35)/5))+ 0.1/(1+exp((V-50)/200))   );
tau_m=@(V)(alpha_m(V)* beta_m(V)   );

%%%%%%%%%%%%%% h %%%%%%%%%%%%%%%%

h_inf=@(V)(   1/(1+exp((V+71.55)/7.43))^2   );
alpha_h=@(V)(   (V<-40)*(0.057*exp(-(V+80)/6.8))   );
beta_h=@(V) (   (V>=-40)*(0.77/(0.13*(1+exp(-(V+10.66)/11.1))))+(V<-40)*(2.7*exp(0.079*V)+3.1*(10^5)*exp(0.3485*V))   );
tau_h=@(V)(   1/(alpha_h(V) + beta_h(V))   );

%%%%%%%%%%%%%% j %%%%%%%%%%%%%%%%

j_inf =@(V)(   1/((1+exp((V+71.55)/7.43))^2)   );
alpha_j= @(V)(   (V<-40)*(-2.5428e4*exp(0.24444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)/(1+exp(0.311*(V+79.23)))   );
beta_j=@(V) (   (V>=-40)*(0.6*exp(0.057*V)/(1+exp(-0.1*(V+32)))) +(V<-40)*(0.02424*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))))  );
tau_j=@(V)(   1/(alpha_j(V) + beta_j(V)));

%%%%%%%%%%%%%%%%%%%%%% d % %%%%%%%%%%%%%%%%%%%

d_inf=@(V) (   1/(1+exp((-V-5)/7.5))   );
alpha_d=@(V) (   (1.4/(1+exp((-V-35)/13)))+0.25   );
beta_d=@(V) (   (1.4/(1+exp((V+5)/5)))   );
gamma_d=@(V) (   (1/(1+exp((50-V)/20)))   );
tau_d=@(V) (   alpha_d(V)*beta_d(V)+gamma_d(V)   );

%%%%%%%%%%%%%%%%%%%%%% f %%%%%%%%%%%%%%%%%%%%

f_inf=@(V) (   1/(1+exp((V+20)/7))   );
%CODICE
tau_f=@(V) (   1125*exp((-(V+27)^2)/300)+(165/(1+exp((25-V)/10)))+80   );
%PAPER
%tau_f=@(V) (   1125*exp((-(V+27)^2)/240)+(165/(1+exp((25-V)/10)))+80   );

%%%%%%%%%%%%%%%%%%%% fca %%%%%%%%%%%%%%%%%%%%

alpha_fca=@(Ca) (   1/(1+(Ca/0.000325)^8)   );
beta_fca=@(Ca) ( 0.1/(1+exp((Ca-0.0005)/0.0001))   );
gamma_fca=@(Ca) ( 0.2/(1+exp((Ca-0.00075)/0.0008))   );
fca_inf=@(Ca) ((alpha_fca(Ca)+beta_fca(Ca)+gamma_fca(Ca)+0.23)/1.46   );
tau_fca= 2 ;

k=@(V,fca, Ca)(  (fca_inf(Ca)<=fca  || V<=threshold)  );    %attenzione quando la integri


%%%%%%%%%%%%%%%%%%%%%% r %%%%%%%%%%%%%%%%%%%

r_inf=@(V) (    1/(1+exp((20-V)/6))   );
tau_r=@(V) (   9.5*exp(-((V+40)^2)/1800)+0.8   );

%%%%%%%%%%%%%%%%%%%%%% s  %%%%%%%%%%%%%%%%%%%


if (epicardium == 1 || Mcell ==1)
    
    s_inf= @(V) (    1/(1+exp((20+V)/5))   );
    tau_s=@(V) (   85*exp(-((V+45)^2)/320)+5/(1+exp((V-20)/5))+3   );
    
elseif ( endo==1 )
    
    s_inf= @(V) (   1/(1+exp((28+V)/5))   );
    tau_s= @(V) (   1000*exp(-(V+67)^2/1000)+8   );
    
end

%%%%%%%%%%%%%%%%%%%%%% chis %%%%%%%%%%%%%%%%%%%

chis_inf=@(V) (   1.0/(1+exp((-5-V)/14))   );
alpha_chis=@(V) (   1100/sqrt(1+exp((-10-V)/6))   );    
beta_chis=@(V) (   1/(1+exp((V-60)/20))   );
tau_chis=@(V) (   alpha_chis(V) * beta_chis(V)    );


%%%%%%%%%%%%%%%%%%%%%% chir1 %%%%%%%%%%%%%%%%%%%

chir1_inf=@(V) (   1/(1+exp((-26.0-V)/7.0))   );
alpha_chir1=@(V) (   450/(1+exp((-45.0-V)/10.0))   );
beta_chir1=@(V) (   6/(1+exp((V+30.0)/11.5))    );
tau_chir1=@(V) (   alpha_chir1(V) * beta_chir1(V)   );


%%%%%%%%%%%%%%%%%%%%%% chir2 %%%%%%%%%%%%%%%%%%%


chir2_inf=@(V) (    1/(1+exp((V+88.0)/24.0))   )   ;
alpha_chir2=@(V) (   3.0/(1+exp((-60-V)/20.0))   )   ;
beta_chir2=@(V) (   1.12/(1+exp((V-60)/20.0))   );
tau_chir2=@(V) (    alpha_chir2(V) * beta_chir2(V)   );



%%%%%%%%%%%%%%%   gating_for_calcium_dynamics  %%%%%

g_inf=@(Ca)(   (Ca<=0.00035)*(1/(1+(Ca/0.00035)^6))+(Ca>0.00035)*(1/(1+(Ca/0.00035)^16))    );

tau_g=2;   %attenzione quando integri

kg=@(V,g, Ca)(  (g_inf(Ca)<=g  || V<=threshold)  );    %attenzione quando la integri


for i=2:length(t)
    
    V_old=V_vect(i-1);
    
    m_old=m_vect(i-1);
    h_old=h_vect(i-1);
    j_old= j_vect(i-1);
    d_old=d_vect(i-1);
    f_old=f_vect(i-1);
    fca_old=fca_vect(i-1);
    r_old=r_vect(i-1);
    s_old=s_vect(i-1);
    chis_old=chis_vect(i-1);
    chir1_old=chir1_vect(i-1);
    chir2_old=chir2_vect(i-1);
    g_old=g_vect(i-1);
    
    Ca_i_old      =Ca_i_vect(i-1);
    Ca_i_tot_old  =Ca_i_tot_vect(i-1);
    Ca_sr_old     =Ca_sr_vect(i-1);
    Ca_sr_tot_old =Ca_sr_tot_vect(i-1);
    Na_i_old      =Na_i_vect(i-1);
    K_i_old       =K_i_vect(i-1);
    
    %%%%%%%%    EE
    
    % variabile m
    L=-1/tau_m(V_old);
    N=m_inf(V_old)/tau_m(V_old);
    [m_new]=integrate_gating(L,N,m_old,tau);
    
    % variabile h
    L=-1/tau_h(V_old);
    N=h_inf(V_old)/tau_h(V_old);
    [h_new]=integrate_gating(L,N,h_old,tau);
    
    % variabile j
    L=-1/tau_j(V_old);
    N=j_inf(V_old)/tau_j(V_old);
    [j_new]=integrate_gating(L,N,j_old,tau);
    
    % variabile d
    L=-1/tau_d(V_old);
    N=d_inf(V_old)/tau_d(V_old);
    [d_new]=integrate_gating(L,N,d_old,tau);
    
    % variabile f
    L=-1/tau_f(V_old);
    N=f_inf(V_old)/tau_f(V_old);
    [f_new]=integrate_gating(L,N,f_old,tau);
    
    % variabile f_Ca
    L=-1/tau_fca;
    N=fca_inf(Ca_i_old)/tau_fca;
    
    L=L*k(V_old,fca_old,Ca_i_old);
    N=N*k(V_old,fca_old,Ca_i_old);
    
    %[fca_new]=integrate_gating(L,N,g_old,tau);
    fca_new = fca_old + tau*(L*fca_old+N);
    
    % variabile r
    L=-1/tau_r(V_old);
    N=r_inf(V_old)/tau_r(V_old);
    [r_new]=integrate_gating(L,N,r_old,tau);
    
    % variabile s
    L=-1/tau_s(V_old);
    N=s_inf(V_old)/tau_s(V_old);
    [s_new]=integrate_gating(L,N,s_old,tau);
    
    % variabile chis
    L=-1/tau_chis(V_old);
    N=chis_inf(V_old)/tau_chis(V_old);
    [chis_new]=integrate_gating(L,N,chis_old,tau);
    
    % variabile chir1
    L=-1/tau_chir1(V_old);
    N=chir1_inf(V_old)/tau_chir1(V_old);
    [chir1_new]=integrate_gating(L,N,chir1_old,tau);
    
    % variabile chir2
    L=-1/tau_chir2(V_old);
    N=chir2_inf(V_old)/tau_chir2(V_old);
    [chir2_new]=integrate_gating(L,N,chir2_old,tau);
    
    % variabile f_Ca
    L=-1/tau_g;
    g_inf(Ca_i_old);
    N=g_inf(Ca_i_old)/tau_g;
    
    L=L*kg(V_old,g_old, Ca_i_old);
    N=N*kg(V_old,g_old, Ca_i_old);
    
    %[g_new]=integrate_gating(L,N,g_old,tau);
    g_new = g_old + tau*(L*g_old+N);
    
    %%%%%%%% COMPUTE CURRENTS OLD   12 + 3
    
    I_Na_old=I_Na(V_old,m_old,h_old,j_old,Na_i_old);
    I_CaL_old=I_CaL(V_old,d_old,f_old,fca_old,Ca_i_old);
    I_to_old=I_to(V_old,r_old,s_old,K_i_old);
    I_Ks_old=I_Ks(V_old,chis_old,Na_i_old,K_i_old);
    I_Kr_old=I_Kr(V_old,chir1_old,chir2_old,K_i_old);
    I_K1_old=I_K1(V_old,K_i_old);
    
    I_NaCa_old=I_NaCa(V_old,Na_i_old,Ca_i_old);
    I_NaK_old=I_NaK(V_old,Na_i_old);
    I_pCa_old=I_pCa(Ca_i_old);
    I_pK_old=I_pK(V_old,K_i_old);
    I_bNa_old=I_bNa(V_old,Na_i_old);
    I_bCa_old=I_bCa(V_old,Ca_i_old);
    
    I_leak_old=I_leak(Ca_i_old,Ca_sr_old);
    I_up_old=I_up(Ca_i_old);
    I_rel_old=I_rel(d_old,g_old,Ca_sr_old);
    
    %%%%%%%%%  INTEGRATE POTENTIAL
    
    I_Ion_old=I_Na_old+I_K1_old+I_to_old+I_Kr_old+I_Ks_old+I_CaL_old+I_NaCa_old+I_NaK_old+I_pCa_old+I_pK_old+I_bNa_old+I_bCa_old;
    V_new=V_old-tau*I_Ion_old; %/Capacitance
   
    
    %%%%%%%%%  INTEGRATE IONIC CONCETRATIONS
    dNa=-(I_Na_old+I_bNa_old+3*I_NaK_old+3*I_NaCa_old)/(V_C *F)*Capacitance;
    Na_i_new= Na_i_old+dNa*tau;

    dK= -(I_K1_old+I_to_old+I_Kr_old+I_Ks_old-2*I_NaK_old+I_pK_old)/(V_C *F)*Capacitance;
    K_i_new=K_i_old+dK*tau;
    
    dCa_sr=(I_up_old-I_rel_old-I_leak_old)*V_C/V_SR;
    Ca_sr_tot_new=Ca_sr_tot_old + tau*dCa_sr;
   
    XT=Ca_sr_tot_new;
    bj=Buf_sr+K_bufsr-XT;
    cj=-XT*K_bufsr;
    
    delta=bj^2-4*cj;
    if (delta<0)
        disp('errore Ca_sr')
        return;
    end
    
    Ca_sr_new= (-bj + sqrt(delta))/2;
    %Ca_sr_new=Ca_sr_old;

    
    dCa_i_tot=-(I_CaL_old+I_bCa_old+I_pCa_old-2*I_NaCa_old)/(2*V_C*F)*Capacitance -I_up_old+I_rel_old+I_leak_old;
    Ca_i_tot_new=Ca_i_tot_old+tau*dCa_i_tot;

    XT=Ca_i_tot_new;
    bj=Buf_c+K_bufc-XT;
    cj=-XT*K_bufc;
    delta=bj^2-4*cj;
    if (delta<0)
        disp('errore Ca_i')
        return;
    end
    
    
    Ca_i_new= (-bj + sqrt(delta))/2;
    %Ca_i_new=Ca_i_old;
    
    
    V_vect(i)=V_new;
    
    m_vect(i)=m_new;
    h_vect(i)=h_new;
    j_vect(i)=j_new;
    chir1_vect(i)=chir1_new; % DIVERSO DA C++ ma dopo il limite di precisione
    chir2_vect(i)=chir2_new;
    chis_vect(i)=chis_new;
    s_vect(i)=s_new;
    r_vect(i)=r_new; % DIVERSO DA C++ ma dopo il limite di precisione
    d_vect(i)=d_new;
    f_vect(i)=f_new;          % DIVERSO DA PAPER
    fca_vect(i)=fca_new;
    g_vect(i)=g_new;          % PROBABILMENTE SBAGLIATO
    
    
    Ca_i_vect(i)=Ca_i_new;
    Ca_i_tot_vect(i)=Ca_i_tot_new ;
    Ca_sr_vect(i)=Ca_sr_new;     
    Ca_sr_tot_vect(i)= Ca_sr_tot_new ;
    Na_i_vect(i)= Na_i_new;
    K_i_vect(i)= K_i_new;
    
end

%plot(t,Ca_i_vect);

end

