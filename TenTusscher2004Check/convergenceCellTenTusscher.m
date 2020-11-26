% oneovertau=[2.^[2:1:8]];
% h=1./oneovertau
% hcube=h.^3

%simulazione predictor corrector
%load('ref_TenTusscher.mat')
%  tic
[Vref,tref, Ca_i_ref, Ca_sr_ref, m_ref,h_ref,j_ref,d_ref,f_ref,r_ref, s_ref,chis_ref,chir1_ref,chir2_ref,g_ref]=CellTenTusscher(2^15,100);
% toc
% 
 [Vcc,tcc, Ca_i_cc, Ca_sr_cc,m_cc,h_cc,j_cc,d_cc,f_cc,r_cc, s_cc,chis_cc,chir1_cc,chir2_cc,g_cc]=CellTenTusscher(2^2,100);
 [Vc,tc, Ca_i_c, Ca_sr_c, m_c,h_c,j_c,d_c,f_c,r_c, s_c,chis_c,chir1_c,chir2_c,g_c]=CellTenTusscher(2^3,100);
 
 [Vm,tm, Ca_i_m, Ca_sr_m,m_m,h_m,j_m,d_m,f_m,r_m, s_m,chis_m,chir1_m,chir2_m,g_m]=CellTenTusscher(2^4,100);
 [Vg,tg, Ca_i_g, Ca_sr_g,m_g,h_g,j_g,d_g,f_g,r_g, s_g,chis_g,chir1_g,chir2_g,g_g]=CellTenTusscher(2^5,100);
 [Vb,tb, Ca_i_b, Ca_sr_b,m_b,h_b,j_b,d_b,f_b,r_b, s_b,chis_b,chir1_b,chir2_b,g_b]=CellTenTusscher(2^6,100);
 [Vr,tr, Ca_i_r, Ca_sr_r,m_r,h_r,j_r,d_r,f_r,r_r, s_r,chis_r,chir1_r,chir2_r,g_r]=CellTenTusscher(2^7,100);
 [Vf,tf, Ca_i_f, Ca_sr_f,m_f,h_f,j_f,d_f,f_f,r_f, s_f,chis_f,chir1_f,chir2_f,g_f]=CellTenTusscher(2^8,100);




% 
%  [Vcc,tcc, Ca_i_cc, Ca_sr_cc, Ca_ss_cc,m_cc,h_cc,j_cc,d_cc,f_cc,f2_cc,fcass_cc,r_cc, s_cc,chis_cc,chir1_cc,chir2_cc,g_cc, f2_inf_cc,fcass_inf_cc,f_inf_cc]=cellTenTusscher2ndOrder(2^2);2
%  [Vc,tc, Ca_i_c, Ca_sr_c, Ca_ss_c,m_c,h_c,j_c,d_c,f_c,f2_c,fcass_c,r_c, s_c,chis_c,chir1_c,chir2_c,g_c, f2_inf_c,fcass_inf_c,f_inf_c]=cellTenTusscher2ndOrder(2^3);3
%  
%  [Vm,tm, Ca_i_m, Ca_sr_m, Ca_ss_m,m_m,h_m,j_m,d_m,f_m,f2_m,fcass_m,r_m, s_m,chis_m,chir1_m,chir2_m,g_m, f2_inf_m,fcass_inf_m,f_inf_m]=cellTenTusscher2ndOrder(2^4);4
%  [Vg,tg, Ca_i_g, Ca_sr_g, Ca_ss_g,m_g,h_g,j_g,d_g,f_g,f2_g,fcass_g,r_g, s_g,chis_g,chir1_g,chir2_g,g_g, f2_inf_g,fcass_inf_g,f_inf_g]=cellTenTusscher2ndOrder(2^5);5
%  [Vb,tb, Ca_i_b, Ca_sr_b, Ca_ss_b,m_b,h_b,j_b,d_b,f_b,f2_b,fcass_b,r_b, s_b,chis_b,chir1_b,chir2_b,g_b, f2_inf_b,fcass_inf_b,f_inf_b]=cellTenTusscher2ndOrder(2^6);6
%  [Vr,tr, Ca_i_r, Ca_sr_r, Ca_ss_r,m_r,h_r,j_r,d_r,f_r,f2_r,fcass_r,r_r, s_r,chis_r,chir1_r,chir2_r,g_r, f2_inf_r,fcass_inf_r,f_inf_r]=cellTenTusscher2ndOrder(2^7);7
%  [Vf,tf, Ca_i_f, Ca_sr_f, Ca_ss_f,m_f,h_f,j_f,d_f,f_f,f2_f,fcass_f,r_f, s_f,chis_f,chir1_f,chir2_f,g_f, f2_inf_f,fcass_inf_f,f_inf_f]=cellTenTusscher2ndOrder(2^8);8
% 




%%%%%%%%%%%%% Action potential %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E1=norm((Vref(1:2^13:end)-Vcc),inf);
E2=norm((Vref(1:2^12:end)-Vc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
E3=norm((Vref(1:2^11:end)-Vm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
E4=norm((Vref(1:2^10:end)-Vg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
E5=norm((Vref(1:2^9:end)-Vb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
E6=norm((Vref(1:2^8:end)-Vr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
E7=norm((Vref(1:2^7:end)-Vf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)


figure(1)
hold on
Errors=[E1,E2,E3,E4,E5,E6,E7]
semilogy(oneovertau(1:7),Errors,'g');
semilogy(oneovertau(1:7),1./oneovertau(1:7).^1,'g');

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

figure(3)
plot(tref(1:50:end),Vref(1:50:end),'k')
hold 
plot(tb,Vb,'b')
plot(tg,Vg,'g')
plot(tm,Vm,'m')
plot(tc,Vc,'c')


% 
% %%%%%%%%%%%%% Ca_i %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% E1=norm((Ca_i_ref(1:2^13:end)-Ca_i_cc),inf);
% E2=norm((Ca_i_ref(1:2^12:end)-Ca_i_c),inf); %to haCa_i_e the absolute add  /norm(Ca_i_ref(1:64:end),inf)
% E3=norm((Ca_i_ref(1:2^11:end)-Ca_i_m),inf); %to haCa_i_e the absolute add /norm(Ca_i_ref(1:32:end),inf)
% E4=norm((Ca_i_ref(1:2^10:end)-Ca_i_g),inf); %to haCa_i_e the absolute add /norm(Ca_i_ref(1:16:end),inf)
% E5=norm((Ca_i_ref(1:2^9:end)-Ca_i_b),inf); %to haCa_i_e the absolute add /norm(Ca_i_ref(1:8:end),inf)
% E6=norm((Ca_i_ref(1:2^8:end)-Ca_i_r),inf); %to haCa_i_e the absolute add /norm(Vref(1:4:end),inf)
% E7=norm((Ca_i_ref(1:2^7:end)-Ca_i_f),inf); %to haCa_i_e the absolute add /norm(Ca_i_ref(1:2:end),inf)
% 
% figure(1)
% hold on
% Errors=[E1,E2,E3,E4,E5,E6,E7]
% loglog(oneovertau(1:7),Errors,'g');
% loglog(oneovertau(1:6),1./oneovertau(1:6).^1,'g');
% 
% %%%% convergence order wrt grid refinement %%%%%%
% 
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% 
% 
% Ratios=[R1,R2,R3,R4,R5,R6]
% figure(2)
% hold on
% semilogy(oneovertau(1:6),Ratios,'g');
% 
% figure(3)
% plot(tref(1:50:end),Ca_i_ref(1:50:end),'k')
% hold 
% plot(tb,Ca_i_b,'b')
% plot(tg,Ca_i_g,'g')
% plot(tm,Ca_i_m,'m')
% plot(tc,Ca_i_c,'c')
% 
% 
% 
% %%%%%%%%%%%%% Ca_sr %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% E1=norm((Ca_sr_ref(1:2^13:end)-Ca_sr_cc),inf);
% E2=norm((Ca_sr_ref(1:2^12:end)-Ca_sr_c),inf); %to haCa_sr_e the absolute add  /norm(Ca_sr_ref(1:64:end),inf)
% E3=norm((Ca_sr_ref(1:2^11:end)-Ca_sr_m),inf); %to haCa_sr_e the absolute add /norm(Ca_sr_ref(1:32:end),inf)
% E4=norm((Ca_sr_ref(1:2^10:end)-Ca_sr_g),inf); %to haCa_sr_e the absolute add /norm(Ca_sr_ref(1:16:end),inf)
% E5=norm((Ca_sr_ref(1:2^9:end)-Ca_sr_b),inf); %to haCa_sr_e the absolute add /norm(Ca_sr_ref(1:8:end),inf)
% E6=norm((Ca_sr_ref(1:2^8:end)-Ca_sr_r),inf); %to haCa_sr_e the absolute add /norm(Vref(1:4:end),inf)
% E7=norm((Ca_sr_ref(1:2^7:end)-Ca_sr_f),inf); %to haCa_sr_e the absolute add /norm(Ca_sr_ref(1:2:end),inf)
% 
% figure(1)
% hold on
% Errors=[E1,E2,E3,E4,E5,E6,E7]
% loglog(oneovertau(1:7),Errors,'g');
% loglog(oneovertau(1:6),1./oneovertau(1:6).^1,'g');
% 
% %%%% convergence order wrt grid refinement %%%%%%
% 
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% 
% 
% Ratios=[R1,R2,R3,R4,R5,R6]
% figure(2)
% hold on
% semilogy(oneovertau(1:6),Ratios,'g');
% 
% figure(3)
% plot(tref(1:50:end),Ca_sr_ref(1:50:end),'k')
% hold 
% plot(tb,Ca_sr_b,'b')
% plot(tg,Ca_sr_g,'g')
% plot(tm,Ca_sr_m,'m')
% plot(tc,Ca_sr_c,'c')
% 
% 
% 
% 
% %%%%%%%%%%%%% Ca_ss %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% E1=norm((Ca_ss_ref(1:2^13:end)-Ca_ss_cc),inf);
% E2=norm((Ca_ss_ref(1:2^12:end)-Ca_ss_c),inf); %to haCa_ss_e the absolute add  /norm(Ca_ss_ref(1:64:end),inf)
% E3=norm((Ca_ss_ref(1:2^11:end)-Ca_ss_m),inf); %to haCa_ss_e the absolute add /norm(Ca_ss_ref(1:32:end),inf)
% E4=norm((Ca_ss_ref(1:2^10:end)-Ca_ss_g),inf); %to haCa_ss_e the absolute add /norm(Ca_ss_ref(1:16:end),inf)
% E5=norm((Ca_ss_ref(1:2^9:end)-Ca_ss_b),inf); %to haCa_ss_e the absolute add /norm(Ca_ss_ref(1:8:end),inf)
% E6=norm((Ca_ss_ref(1:2^8:end)-Ca_ss_r),inf); %to haCa_ss_e the absolute add /norm(Vref(1:4:end),inf)
% E7=norm((Ca_ss_ref(1:2^7:end)-Ca_ss_f),inf); %to haCa_ss_e the absolute add /norm(Ca_ss_ref(1:2:end),inf)
% 
% figure(1)
% hold on
% Errors=[E1,E2,E3,E4,E5,E6,E7]
% loglog(oneovertau(1:7),Errors,'g');
% loglog(oneovertau(1:6),1./oneovertau(1:6).^1,'g');
% 
% %%%% convergence order wrt grid refinement %%%%%%
% 
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% 
% 
% Ratios=[R1,R2,R3,R4,R5,R6]
% figure(2)
% hold on
% semilogy(oneovertau(1:6),Ratios,'g');
% 
% figure(3)
% plot(tref(1:50:end),Ca_ss_ref(1:50:end),'k')
% hold 
% plot(tb,Ca_ss_b,'b')
% plot(tg,Ca_ss_g,'g')
% plot(tm,Ca_ss_m,'m')
% %plot(tc,Ca_ss_c,'c') %oscilla


%%%%%%%%%%%%% m gate %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E1=norm((m_ref(1:2^13:end)-m_cc),inf);
E2=norm((m_ref(1:2^12:end)-m_c),inf); %to ham_e the absolute add  /norm(m_ref(1:64:end),inf)
E3=norm((m_ref(1:2^11:end)-m_m),inf); %to ham_e the absolute add /norm(m_ref(1:32:end),inf)
E4=norm((m_ref(1:2^10:end)-m_g),inf); %to ham_e the absolute add /norm(m_ref(1:16:end),inf)
E5=norm((m_ref(1:2^9:end)-m_b),inf); %to ham_e the absolute add /norm(m_ref(1:8:end),inf)
E6=norm((m_ref(1:2^8:end)-m_r),inf); %to ham_e the absolute add /norm(Vref(1:4:end),inf)
E7=norm((m_ref(1:2^7:end)-m_f),inf); %to ham_e the absolute add /norm(m_ref(1:2:end),inf)

figure(1)
hold on
Errors=[E1,E2,E3,E4,E5,E6,E7]
loglog(oneovertau(1:7),Errors,'g');
loglog(oneovertau(1:6),1./oneovertau(1:6).^1,'g');

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

figure(3)
plot(tref(1:50:end),m_ref(1:50:end),'k')
hold 
plot(tb,m_b,'b')
plot(tg,m_g,'g')
plot(tm,m_m,'m')
plot(tc,m_c,'c')
% 
% %%%%%%%%%%%%% f gate %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% E1=norm((fcass_ref(1:2^13:end)-fcass_cc),inf);
% E2=norm((fcass_ref(1:2^12:end)-fcass_c),inf); %to hafcass_e the absolute add  /norm(fcass_ref(1:64:end),inf)
% E3=norm((fcass_ref(1:2^11:end)-fcass_m),inf); %to hafcass_e the absolute add /norm(fcass_ref(1:32:end),inf)
% E4=norm((fcass_ref(1:2^10:end)-fcass_g),inf); %to hafcass_e the absolute add /norm(fcass_ref(1:16:end),inf)
% E5=norm((fcass_ref(1:2^9:end)-fcass_b),inf); %to hafcass_e the absolute add /norm(fcass_ref(1:8:end),inf)
% E6=norm((fcass_ref(1:2^8:end)-fcass_r),inf); %to hafcass_e the absolute add /norm(Vref(1:4:end),inf)
% E7=norm((fcass_ref(1:2^7:end)-fcass_f),inf); %to hafcass_e the absolute add /norm(fcass_ref(1:2:end),inf)
% 
% figure(1)
% hold on
% Errors=[E1,E2,E3,E4,E5,E6,E7]
% loglog(oneovertau(1:7),Errors,'g');
% loglog(oneovertau(1:6),1./oneovertau(1:6).^1,'g');
% 
% %%%% convergence order wrt grid refinement %%%%%%
% 
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% 
% 
% Ratios=[R1,R2,R3,R4,R5,R6]
% figure(2)
% hold on
% semilogy(oneovertau(1:6),Ratios,'g');
% 
% figure(3)
% plot(tref(1:50:end),fcass_ref(1:50:end),'k')
% hold 
% plot(tb,fcass_b,'b')
% plot(tg,fcass_g,'g')
% plot(tm,fcass_m,'m')
% plot(tc,fcass_c,'c')

