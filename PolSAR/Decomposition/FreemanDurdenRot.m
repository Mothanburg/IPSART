% FreemanDurdenRot - 去定向的Freeman-Durden三分量分解
% 参考：10.1109/TGRS.2010.2041242
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART 
function [Ps,Pd,Pv] = FreemanDurdenRot(T3)

arguments
    T3 PolT3
end

height = T3.Height;
width = T3.Width;

B = (T3.m22 - T3.m33) / 2;
E = T3.m23_r; % real(T3.PageAt(2, 3) + T3.PageAt(3, 2)) / 2;
cos4t = B ./ sqrt(B.^2 + E.^2);
sin4t = E ./ sqrt(B.^2 + E.^2);
cos2t = sqrt((1 + cos4t) / 2);
sin2t = sin4t ./ (2 * cos2t);

Ps = zeros(height, width, T3.Dtype);
Pd = zeros(height, width, T3.Dtype);
Pv = zeros(height, width, T3.Dtype);

T3 = parallel.pool.Constant(T3);
parfor j = 1:width
    for i = 1:height
        q = [1 0 0; 0 cos2t(i,j) sin2t(i,j); 0 -sin2t(i,j) cos2t(i,j)];
        t = T3.Value.MatAt(i, j);
        t1 = q * t * q';
        if t1(1,1) <= t1(3,3)
            Pv(i,j) = real(3 * t1(1,1));
            Ps(i,j) = 0;
            Pd(i,j) = real(t1(2,2) + t1(3,3) - 2 * t1(1,1));
        else
            Pv(i,j) = 3 * real(t1(3,3));
            x11 = real(t1(1,1) - t1(3,3));
            x22 = real(t1(2,2) - t1(3,3));
            if abs(t1(1,2))^2 > x11 * x22
                if x11 > x22
                    Ps(i,j) = x11 + x22;
                    Pd(i,j) = 0;
                else
                    Ps(i,j) = 0;
                    Pd(i,j) = x11 + x22;
                end
            else
                if x11 > x22
                    Ps(i,j) = x11 + abs(t1(1,2))^2 / x11;
                    Pd(i,j) = x22 - abs(t1(1,2))^2 / x11;
                else
                    Ps(i,j) = x11 - abs(t1(1,2))^2 / x22;
                    Pd(i,j) = x22 + abs(t1(1,2))^2 / x22;
                end
            end
        end
    end
end

end
