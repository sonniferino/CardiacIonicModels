%load('../fine_data_AB2AM3_also_gating.mat'); 
oppure 

[Vref,tref,vref,mref,fref,toref,Xref]=bernus_initialized_gating(num_steps)


%tentativo di interpolazione 
% tused=tref([1,828,1494,9526:10000:end,end])
% Vused=Vref([1,828,1494,9526:10000:end,end])
% p=polyfit(tused,Vused,5);  %p(1)x^4+.....+p(5)
% 
% 
% V=polyval(p,tref);
% 
% close all
% plot(tref,V,tref,Vref)



