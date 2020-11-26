oneovertau=[2.^[2:1:8]];
h=1./oneovertau
hcube=h.^3

%simulazione predictor corrector

load('fine_data_AB2AM3_also_gating.mat') 
 tic
%[Vref,tref, I_Naref,I_Caref,I_toref, I_Kref,I_K1ref,I_Nabref,I_Cabref,I_NaKref,I_NaCaref]=bernus_AB2AM3_also_gating(2^15); 
 toc

 [Vcc,tcc]=bernus_RL(2^2);
 [Vc,tc]=bernus_RL(2^3);
 
 [Vm,tm]=bernus_RL(2^4);
 [Vg,tg]=bernus_RL(2^5);
 [Vb,tb]=bernus_RL(2^6);
 [Vr,tr]=bernus_RL(2^7);
 [Vf,tf]=bernus_RL(2^8);


%%%%% alternative extrapolate and use L infty
E1=norm((Vref(1:2^13:end)-Vcc),inf);
E2=norm((Vref(1:2^12:end)-Vc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
E3=norm((Vref(1:2^11:end)-Vm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
E4=norm((Vref(1:2^10:end)-Vg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
E5=norm((Vref(1:2^9:end)-Vb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
E6=norm((Vref(1:2^8:end)-Vr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
E7=norm((Vref(1:2^7:end)-Vf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)

%%%%%% alternative exrtapolate and use relative root square mean error RRMS

% 
% E1=sqrt((1/length(Vc))*sum(((Vref(1:64:end)-Vc).^2))/sum(((Vref(1:64:end)).^2)));
% E2=sqrt((1/length(Vm))*sum(((Vref(1:32:end)-Vm).^2))/sum(((Vref(1:32:end)).^2)));
% E3=sqrt((1/length(Vg))*sum(((Vref(1:16:end)-Vg).^2))/sum(((Vref(1:16:end)).^2)));
% E4=sqrt((1/length(Vb))*sum(((Vref(1:8:end)-Vb).^2))/sum(((Vref(1:8:end)).^2)));
% E5=sqrt((1/length(Vr))*sum(((Vref(1:4:end)-Vr).^2))/sum(((Vref(1:4:end)).^2)));
% E6=sqrt((1/length(Vf))*sum(((Vref(1:2:end)-Vf).^2))/sum(((Vref(1:2:end)).^2)));

%%%%%% plot the error %%%%%

figure(1)
hold on
Errors=[E1,E2,E3,E4,E5,E6,E7]
semilogy(oneovertau(1:7),Errors,'g');

%%%% convergence order wrt grid refinement %%%%%%

R1=log2(E1/E2);
R2=log2(E2/E3);
R3=log2(E3/E4);
R4=log2(E4/E5);
R5=log2(E5/E6);
R6=log2(E6/E7);


Ratios=[R1,R2,R3,R4,R5,R6]
figure(2)
hold on
semilogy(oneovertau(1:6),Ratios,'g');
