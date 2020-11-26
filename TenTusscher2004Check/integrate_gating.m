function [s_new]=integrate_gating(L,N,s_old,dt)

         %if (abs(L)<1e-30 && abs(N)<1e-30)
         %    s_new=s_old;
         %else
             s_new   = exp(L*dt)*s_old+1/L*(exp(L*dt)-1)*N; 
         %end
end