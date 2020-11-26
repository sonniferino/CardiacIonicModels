Capacitance=0.185; %nel paper

K_o=5.4; %controllato paper e codice
Na_o=140.0; %nel paper
Ca_o= 2.0;%controllato: loro codice

G_Na=14.838;   %controllato: loro codice e paper
G_CaL= 0.000175;    %controllato: loro codice nel paper e' invece (1.75)^(-4);

if (epicardium==1)
G_to=0.294;  %controllato: loro codice e paper
else
G_to=0.073; %controllato: loro codice e codice
end


G_Ks=0.245;  % fine for epi/endo cells %controllato
G_Kr= 0.096; %controllato paper e codice 
G_K1= 5.405; %controllato paper e codice

gamma= 0.35; %corretto paper e codice
alpha= 2.5; %controllato: paper e codice

k_NaCa= 1000;
K_mNai= 87.5; %controllato: paper e  codice
K_mCa= 1.38;%controllato: paper e codice
k_sat= 0.1;%controllato: paper e codice

P_NaK= 1.362; %corretto: paper e codice

K_mK= 1.0;   %corretto: paper e codice
K_mNa= 40.0; %corretto: paper e codice
%G_pCa= 0.025; %coerente con paper 0.025
G_pCa= 0.825; %coerente con codice
K_pCa=0.0005;%corretto: paper e codice
G_pK=0.0146;%corretto: paper e codice

G_bNa=0.00029; %corretto: paper e codice
G_bCa=0.000592;%corretto: paper e codice


%V_C=16404; %nel paper 
%V_SR=1094.0; %nel paper
V_C=0.016404; % CODICE
V_SR=0.001094;% CODICE

V_maxup=0.000425; %corretto codice e anche paper
K_up=0.00025; %corretto codice e anche paper
V_leak=0.00008; %corretto

%a_rel=16.464;     % paper sono diverse le unit? di misura
a_rel = 0.016464;  % codice con unit? corrette

b_rel = 0.25;      % valore paper. Nel codice ? gi? messo al quadrato

%c_rel=8.232;      %  paper sono diverse le unit? di misura
c_rel = 0.008232;  % codice con unit? corrette


Buf_sr=10.0;  %corretto paper e codice
K_bufsr=0.3;  %corretto paper e codice

Buf_c= 0.15; %corretto paper e codice
K_bufc= 0.001; % corretto paper e codice
