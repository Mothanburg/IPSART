% Model-based decomposition for dual-pol
% 10.1109/TGRS.2021.3137588
function [ms,mv,alpha,delta] = DualPolModelBased(C2, polTx)

arguments
    C2 PolC2
    polTx string % The transfer polarization, must be "H" or "V"
end

switch upper(polTx)
    case "H"
        polTx = "H";
    case "V"
        polTx = "V";
    otherwise
        error("Invalid transfer polarization '%s'", polTx);
end

height = C2.Height;
width = C2.Width;

s1 = C2.SPAN;
s2 = C2.m11 - C2.m22;
s3 = 2 * C2.m12_r;
s4 = 2 * C2.m12_i;

a = 0.75;
if polTx == "H"
    b = -2 * s1 + 0.5 * s2;
else
    b = -2 * s1 - 0.5 * s2;
end
c = s1.^2 - s2.^2 - s3.^2 - s4.^2;

% "mv" satisfy "a*mv^2 + b*mv + c = 0"
% And "mv" must be greater than 0 and less than "s1"
root1 = (-b + sqrt(b.^2 - 4 * a .* c)) ./ (2 * a);
root2 = (-b - sqrt(b.^2 - 4 * a .* c)) ./ (2 * a);

mv = zeros(height, width, C2.Dtype);
cond1 = root1 < s1;
cond2 = ~cond1 & (root2 >= 0);
mv(cond1) = root1(cond1);
mv(cond2) = root2(cond2);


ms = s1 - mv;
if polTx == "H"
    s2p = s2 - 0.5 * mv;
else
    s2p = s2 + 0.5 * mv;
end
s3p = s3 ./ ms;
s4p = s4 ./ ms;

alpha = acosd(s2p ./ ms) / 2;
alpha(isnan(alpha)) = 0;
alpha = real(alpha);

delta = rad2deg(angle(s3p + 1i * s4p));

end
