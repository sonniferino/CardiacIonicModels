oneovertau=[2.^[0:1:11]];
h=1./oneovertau
hcube=h.^3

% %simulazione predictor corrector
% 
%load('../fine_data_AB2AM3_also_gating.mat') 
%  tic
% %[Vref,tref, I_Naref,I_Caref,I_toref, I_Kref,I_K1ref,I_Nabref,I_Cabref,I_NaKref,I_NaCaref]=bernus_AB2AM3_also_gating(2^15); 
%  toc
%
[Vi,ti,vi,mi,fi,toi,Xi,Lv,Lm,Lf,Lto,LX,Nv,Nm,Nf,Nto,NX]=bernus_IFAB2(2^0,tref,Vref);
% [Vccc,tccc,vccc,mccc,fccc,toccc,Xccc]=bernus_IFAB2(2^1,tref,Vref);
% [Vcc,tcc,vcc,mcc,fcc,tocc,Xcc]=bernus_IFAB2(2^2,tref,Vref);
% [Vc,tc,vc,mc,fc,toc,Xc]=bernus_IFAB2(2^3,tref,Vref);
% [Vm,tm,vm,mm,fm,tom,Xm]=bernus_IFAB2(2^4,tref,Vref);
% [Vg,tg,vg,mg,fg,tog,Xg]=bernus_IFAB2(2^5,tref,Vref);
% [Vb,tb,vb,mb,fb,tob,Xb]=bernus_IFAB2(2^6,tref,Vref);
% [Vr,tr,vr,mr,fr,tor,Xr]=bernus_IFAB2(2^7,tref,Vref);
% [Vf,tf,vf,mf,ff,tof,Xf]=bernus_IFAB2(2^8,tref,Vref);
 [Vff,tff,vff,mff,fff,toff,Xff]=bernus_IFAB2(2^9,tref,Vref);
% [Vfff,tfff,vfff,mfff,ffff,tofff,Xfff]=bernus_IFAB2(2^10,tref,Vref);
% [Vsf,tsf,vsf,msf,fsf,tosf,Xsf]=bernus_IFAB2(2^11,tref,Vref);
% 

% figure(2)
%  Vsupp=Vref(1:20:end);
% 
% subplot(2,5,1)
% plot(Vsupp,Lv)
% legend('L_v')
% 
% subplot(2,5,2)
% plot(Vsupp,Lm)
% legend('L_m')
% 
% subplot(2,5,3)
% plot(Vsupp,Lf)
% legend('L_f')
% subplot(2,5,4)
% plot(Vsupp,Lto)
% legend('L_{to}')
% subplot(2,5,5)
% plot(Vsupp,LX)
% legend('L_x')
% 
% 
% subplot(2,5,6)
% plot(Vsupp,Nv)
% legend('N_v')
% 
% subplot(2,5,7)
% plot(Vsupp,Nm)
% legend('N_m')
% 
% subplot(2,5,8)
% plot(Vsupp,Nf)
% legend('N_f')
% subplot(2,5,9)
% plot(Vsupp,Nto)
% legend('N_{to}')
% subplot(2,5,10)
% plot(Vsupp,NX)
% legend('N_x')


figure(2)
 Vsupp=Vref(1:20:end);
 tsupp=tref(1:20:end);

subplot(1,2,1)
plot(tsupp,Vsupp)
legend('V_ref')

subplot(1,2,2)
plot(tff,vff)
hold on 
plot(tff,mff)
plot(tff,fff)
plot(tff,toff)
plot(tff,Xff)


% 
% 
% %%%%%%%%%%%%%%%%%%%%%% per v %%%%%%%%%%%%%%%%
%  
% %%%%% alternative extrapolate and use L infty
% E0=norm((vref(1:2^15:end)-vi),inf);
% E1=norm((vref(1:2^14:end)-vccc),inf);
% E2=norm((vref(1:2^13:end)-vcc),inf);
% E3=norm((vref(1:2^12:end)-vc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
% E4=norm((vref(1:2^11:end)-vm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
% E5=norm((vref(1:2^10:end)-vg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
% E6=norm((vref(1:2^9:end)-vb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
% E7=norm((vref(1:2^8:end)-vr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E8=norm((vref(1:2^7:end)-vf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E9=norm((vref(1:2^6:end)-vff),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E10=norm((vref(1:2^5:end)-vfff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E11=norm((vref(1:2^4:end)-vsf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% 
% 
% %%%%%% alternative exrtapolate and use relative root square mean error RRMS
% %%%%%% plot the error %%%%%
% figure(1)
% subplot(2,5,1)
% hold on
% Errors=[E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11]
% semilogy(oneovertau(1:12),Errors,'r');
% %%legend('Error for gating v')
% 
% %%%% convergence order wrt grid refinement %%%%%%
% R0=log2(E0/E1);
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% R7=log2(E7/E8);
% R8=log2(E8/E9);
% R9=log2(E9/E10);
% R10=log2(E10/E11);
% 
% 
% subplot(2,5,6)
% hold on
% Ratios=[R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10]
% semilogy(oneovertau(1:11),Ratios,'r');
% 
% %legend('Order of convergence for gating v')
% 
%  %%%%%%%%%%%%%%%%%%%%%%% per m %%%%%%%%%%%%%%%%
%  
% %%%% alternative extrapolate and use L infty
% E0=norm((mref(1:2^15:end)-mi),inf);
% E1=norm((mref(1:2^14:end)-mccc),inf);
% E2=norm((mref(1:2^13:end)-mcc),inf);
% E3=norm((mref(1:2^12:end)-mc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
% E4=norm((mref(1:2^11:end)-mm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
% E5=norm((mref(1:2^10:end)-mg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
% E6=norm((mref(1:2^9:end)-mb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
% E7=norm((mref(1:2^8:end)-mr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E8=norm((mref(1:2^7:end)-mf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E9=norm((mref(1:2^6:end)-mff),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E10=norm((mref(1:2^5:end)-mfff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E11=norm((mref(1:2^4:end)-msf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% 
% %%%%% alternative exrtapolate and use relative root square mean error RRMS
% %%%%% plot the error %%%%%
% figure(1)
% subplot(2,5,2)
% hold on
% Errors=[E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11]
% semilogy(oneovertau(1:12),Errors,'r');
% %legend('Error for gating v')
% 
% %%% convergence order wrt grid refinement %%%%%%
% R0=log2(E0/E1);
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% R7=log2(E7/E8);
% R8=log2(E8/E9);
% R9=log2(E9/E10);
% R10=log2(E10/E11);
% 
% subplot(2,5,7)
% hold on
% Ratios=[R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10]
% semilogy(oneovertau(1:11),Ratios,'r');
% 
% legend('Order of convergence for gating v')
% 
%  %%%%%%%%%%%%%%%%%%%%%%% per f %%%%%%%%%%%%%%%%
%  
% %%%%% alternative extrapolate and use L infty
% E0=norm((fref(1:2^15:end)-fi),inf);
% E1=norm((fref(1:2^14:end)-fccc),inf);
% E2=norm((fref(1:2^13:end)-fcc),inf);
% E3=norm((fref(1:2^12:end)-fc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
% E4=norm((fref(1:2^11:end)-fm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
% E5=norm((fref(1:2^10:end)-fg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
% E6=norm((fref(1:2^9:end)-fb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
% E7=norm((fref(1:2^8:end)-fr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E8=norm((fref(1:2^7:end)-ff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E9=norm((fref(1:2^6:end)-fff),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E10=norm((fref(1:2^5:end)-ffff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E11=norm((fref(1:2^4:end)-fsf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% 
% %%%%%% alternative exrtapolate and use relative root square mean error RRMS
% %%%%%% plot the error %%%%%
% figure(1)
% subplot(2,5,3)
% hold on
% Errors=[E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11]
% semilogy(oneovertau(1:12),Errors,'r');
% %%legend('Error for gating v')
% 
% %%%% convergence order wrt grid refinement %%%%%%
% R0=log2(E0/E1);
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% R7=log2(E7/E8);
% R8=log2(E8/E9);
% R9=log2(E9/E10);
% R10=log2(E10/E11);
% 
% subplot(2,5,8)
% hold on
% Ratios=[R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10]
% semilogy(oneovertau(1:11),Ratios,'r');
% 
% %legend('Order of convergence for gating v')
% 
%  
%  %%%%%%%%%%%%%%%%%%%%%%% per to %%%%%%%%%%%%%%%%
%  
% %%%%% alternative extrapolate and use L infty
% E0=norm((toref(1:2^15:end)-toi),inf);
% E1=norm((toref(1:2^14:end)-toccc),inf);
% E2=norm((toref(1:2^13:end)-tocc),inf);
% E3=norm((toref(1:2^12:end)-toc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
% E4=norm((toref(1:2^11:end)-tom),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
% E5=norm((toref(1:2^10:end)-tog),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
% E6=norm((toref(1:2^9:end)-tob),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
% E7=norm((toref(1:2^8:end)-tor),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E8=norm((toref(1:2^7:end)-tof),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E9=norm((toref(1:2^6:end)-toff),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E10=norm((toref(1:2^5:end)-tofff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E11=norm((toref(1:2^4:end)-tosf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% 
% %%%%%% alternative exrtapolate and use relative root square mean error RRMS
% %%%%%% plot the error %%%%%
% figure(1)
% subplot(2,5,4)
% hold on
% Errors=[E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11]
% semilogy(oneovertau(1:12),Errors,'r');
% %%legend('Error for gating v')
% 
% %%%% convergence order wrt grid refinement %%%%%%
% R0=log2(E0/E1);
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% R7=log2(E7/E8);
% R8=log2(E8/E9);
% R9=log2(E9/E10);
% R10=log2(E10/E11);
% 
% subplot(2,5,9)
% hold on
% Ratios=[R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10]
% semilogy(oneovertau(1:11),Ratios,'r');
% 
% %legend('Order of convergence for gating v')
% 
%  %%%%%%%%%%%%%%%%%%%%%%% per X %%%%%%%%%%%%%%%%
%  
% %%%%% alternative extrapolate and use L infty
% E0=norm((Xref(1:2^15:end)-Xi),inf);
% E1=norm((Xref(1:2^14:end)-Xccc),inf);
% E2=norm((Xref(1:2^13:end)-Xcc),inf);
% E3=norm((Xref(1:2^12:end)-Xc),inf); %to have the absolute add  /norm(Vref(1:64:end),inf)
% E4=norm((Xref(1:2^11:end)-Xm),inf); %to have the absolute add /norm(Vref(1:32:end),inf)
% E5=norm((Xref(1:2^10:end)-Xg),inf); %to have the absolute add /norm(Vref(1:16:end),inf)
% E6=norm((Xref(1:2^9:end)-Xb),inf); %to have the absolute add /norm(Vref(1:8:end),inf)
% E7=norm((Xref(1:2^8:end)-Xr),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E8=norm((Xref(1:2^7:end)-Xf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E9=norm((Xref(1:2^6:end)-Xff),inf); %to have the absolute add /norm(Vref(1:4:end),inf)
% E10=norm((Xref(1:2^5:end)-Xfff),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% E11=norm((Xref(1:2^4:end)-Xsf),inf); %to have the absolute add /norm(Vref(1:2:end),inf)
% 
% %%%%%% alternative exrtapolate and use relative root square mean error RRMS
% %%%%%% plot the error %%%%%
% figure(1)
% subplot(2,5,5)
% hold on
% Errors=[E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11]
% semilogy(oneovertau(1:12),Errors,'r');
% %%legend('Error for gating v')
% 
% %%%% convergence order wrt grid refinement %%%%%%
% R0=log2(E0/E1);
% R1=log2(E1/E2);
% R2=log2(E2/E3);
% R3=log2(E3/E4);
% R4=log2(E4/E5);
% R5=log2(E5/E6);
% R6=log2(E6/E7);
% R7=log2(E7/E8);
% R8=log2(E8/E9);
% R9=log2(E9/E10);
% R10=log2(E10/E11);
% 
% 
% subplot(2,5,10)
% hold on
% Ratios=[R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10]
% semilogy(oneovertau(1:11),Ratios,'r');
% 
% %legend('Order of convergence for gating v')