function [Ps,Pd,Pv] = FreemanDurden4CTLR(G)

arguments
    G (:,:,4)
end

G = real(G);
g0 = G(:,:,1);
g1 = G(:,:,2);
g2 = G(:,:,3);
g3 = G(:,:,4);
x1 = g0-sqrt(g1.^2+g2.^2+g3.^2);
x = x1*0.65;

mask1 = g3<0;
pd1 = ((g0+g3-x).*(g0-g3-x)-g1.^2-g2.^2)./(2*(g0-g3+x));
ps1 = ((g0-g3-x).^2+g1.^2+g2.^2)./(2*(g0-g3-x));
mask2 = ~mask1;
pd2 = ((g0+g3-x).^2+g1.^2+g2.^2)./(2*(g0+g3-x));
ps2 = ((g0-g3-x).*(g0+g3-x)-g1.^2-g2.^2)./(2*(g0+g3-x));

Pv = x;
Pd = pd1.*mask1+pd2.*mask2;
Ps = ps1.*mask1+ps2.*mask2;

end