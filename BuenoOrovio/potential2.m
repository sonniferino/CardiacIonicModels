%clear all
%close all

t=linspace(0,300,19201);

tau=t(2)-t(1);

u=zeros(size(t));
v=zeros(size(t));
w=zeros(size(t));
s=zeros(size(t));

counter=1;

u(counter)=0.5;
v(counter)=1;
w(counter)=1;
s(counter)=0;

clear counter

for i=2:length(t)
    
    uold=u(i-1);
    vold=v(i-1);
    wold=w(i-1);
    sold=s(i-1);
    
    [vnew,wnew,snew]=differentialStep(vold,wold,sold,uold,tau)

    v(i)=vnew;
    w(i)=wnew;
    s(i)=snew;
    
    [J_fi, J_so, J_si]=currents(vnew,wnew,snew,uold);
    
    unew=uold-tau*(J_fi+  J_so + J_si)
    
    u(i)=unew;
    
end

figure(1)
hold on
plot(t,u,'b');%,t,s,t,w,t,v)