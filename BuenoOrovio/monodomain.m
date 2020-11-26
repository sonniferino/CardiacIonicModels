clear all
close all

% aviOutput=fullfile('Movie_bueno_orovio.avi');
% aviobj = VideoWriter(aviOutput);
% aviobj.Quality= 100;
% aviobj.FrameRate = 24;
% open(aviobj);

t=linspace(0,1000,10000);
tau=t(2)-t(1);

x=linspace(0,2,1000);
u=zeros(size(x))';
u(1:50)=1;
v=ones(size(x))';
w=ones(size(x))';
s=zeros(size(x))';
Iion=zeros(size(x))';


h=diff(x)';
M=1/3*diag([h;0]+[0;h])+1/6*diag(h,1)+1/6*diag(h,-1);
A=diag([1./h;0]+[0;1./h])-diag(1./h,1)-diag(1./h,-1);

for i=2:6000
    
    for j=1:length(x)
        
    uold=u(j);
    vold=v(j);
    wold=w(j);
    sold=s(j);
    
    [vold,wold,sold]=differentialStep(vold,wold,sold,uold,tau);

    v(j)=vold;
    w(j)=wold;
    s(j)=sold;
    
    [J_fi, J_so, J_si]=currents(vold,wold,sold,uold);
    
    Iion(j)=-(J_fi+  J_so + J_si);
    
    end

    u=(M+tau*0.01*A)\(M*u+tau*M*Iion);
   
    plot(u)
    ylim([-0.1 1.5])
    title(num2str(i*tau))
 %   %pause(0.01)

%    f=getframe(gcf);
  %  writeVideo(aviobj,f);
     
end


%close(aviobj);

