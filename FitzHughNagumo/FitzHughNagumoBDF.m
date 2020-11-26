


clear all
close all
 
t=linspace(0,1000,100000);
 
V(1)=-80;
w(1)=0.0;
 
counter=1;
 
tau=t(2)-t(1);
Iapp=0;

integration_order=1
 

a = 0.13;
mu1 = 0.26;
mu2 = 0.1 ;
mu3 = 0.013;
mu4 = 1;
 
f=@(t)(t<=0.2);
 
for i=2:length(t)
        Vold=V(i-1);
        wold=w(i-1);

         Iapp=0.000000350*f(t(i-1));

        Vold_norm=(Vold+80)/110;
        
        wnew=wold+ tau*mu3*(Vold_norm-mu4*wold);

        Vnew=Vold+tau*(Iapp-mu1*(Vold-30)*(Vold+80)*(Vold+50) -mu2*wold*(Vold_norm));

        %Vnew_norm=(Vnew+80)/110;
        
        
        
        
         caso = min(counter, integration_order);

        
%         if (caso==1)
% 
%         wnew = (wold +  tau * mu3*( (Vnew_norm)))/(1.0 + tau* mu3 * mu4 );
% 
% 
%         elseif (caso==2)
% 
% 
%         wnew = ( 4.0/3.0 * w(i-1)- 1.0/3.0 * w(i-2) + 2.0/3.0 *  ( tau * mu3*( (Vnew_norm))))/( 1.0 +2.0/3.0 *tau* mu3 * mu4 );
% 
%         elseif (caso==3)
% 
%         wnew = ( 18.0/11.0 * w(i-1)-9.0/11.0 * w(i-2) + 2.0/11.0 * w(i-3)+ 6.0/11.0 * ( tau * mu3*( (Vnew_norm))))/( 1.0 + 6.0/11.0 *tau* mu3 * mu4 );
% 
% 
%         elseif (caso==4)
% 
%         wnew = ( 48.0/25.0 * w(i-1)-36.0/25.0 * w(i-2) + 16.0/25.0 * w(i-3)-3.0/25.0 * w(i-4)+ 12.0/25.0 * ( tau * mu3*( (Vnew_norm))))/( 1.0 + 12.0 / 25.0 *tau* mu3 * mu4 );
% 
%         end 
%         
        V(i)=Vnew;
        w(i)=wnew;
        counter=counter+1;

    
end
 

hold on 
plot(t,V);

figure(2)
plot(t,w);

