close all

run=1;

if (run==1)
clear all

init=5;
fin=14;

for i=init+1:fin
    i
    tic
    [V,t]=CellTenTusscher(2^i,300);
    VV{i-init}=V;
    tt{i-init}=t;
    toc
end

end


for i=1:length(VV)-1
       ratio= ( length( VV{length(VV)} )-1 )./( length( VV{i} )-1 );
       
       err(i)=max(abs( VV{end}(1:ratio:end) - VV{i} ));
       
end

log2(err(1:end-1)./err(2:end))

figure
hold on

for i=1:length(VV)
    
    plot(tt{i},VV{i})
end
