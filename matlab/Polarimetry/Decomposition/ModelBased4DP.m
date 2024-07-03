function [ms,mv,alpha,delta] = ModelBased4DP(C2, polTx)

arguments
    C2 PolC2
    polTx string
end

[ms,mv,alpha,delta] = MBDP_matlab(C2, polTx);

end

function [ms,mv,alpha,delta] = MBDP_matlab(C2, polTx)

height = C2.Height;
width = C2.Width;

s1 = C2.SPAN;
s2 = C2.getPageAt(1, 1) - C2.getPageAt(2, 2);
s3 = 2 * real(C2.getPageAt(1, 2));
s4 = 2 * imag(C2.getPageAt(1, 2));

a = 0.75;
if polTx == "H"
    b = -2 * s1 + 0.5 * s2;
else
    b = -2 * s1 - 0.5 * s2;
end
c = s1.^2 - s2.^2 - s3.^2 - s4.^2;

x1 = (-b + sqrt(b.^2 - 4 * a .* c)) ./ (2 * a);
x2 = (-b - sqrt(b.^2 - 4 * a .* c)) ./ (2 * a);

mv = zeros(height, width, C2.Dtype);
parfor j=1:width
    for i=1:height
        if x1(i,j) < s1(i,j)
            mv(i,j) = x1(i,j);
        else
            if x2(i,j) < 0
                mv(i,j) = 0;
            else
                mv(i,j) = x2(i,j);
            end
        end
    end
end

ms = s1 - mv;
if polTx == "H"
    s2p = (s2 - 0.5 * mv) ./ ms;
else
    s2p = (s2 + 0.5 * mv) ./ ms;
end
s3p = s3 ./ ms;
s4p = s4 ./ ms;

alpha = acosd(-s2p) / 2;
alpha(isnan(alpha)) = 0;
alpha = real(alpha);
delta = rad2deg(angle(s3p + 1i * s4p));

end
