function [FAC_NH] = FAC_NH_func(b,T,srho)
ps = @(m,rho) max(0,m.*(1./rho -1/873));

FAC_NH = NaN(size(b));
avg_rho = NaN(size(b));

for i = 1:size(b,1)
    for j = 1:size(b,2)
        thickness = 0.1;
        
        [rho,~, ~] =  densitymodel(T(i,j),b(i,j),srho,0:thickness:10,...
        'NabarroHerring');
        z_weq = [0 cumsum(rho(1:end-1)*thickness)];

        thickness_weq = z_weq;
        thickness_weq(2:end) = z_weq(2:end) - z_weq(1:end-1);
        FAC_NH(i,j)  = sum(ps(thickness_weq,rho));
    end
end
end