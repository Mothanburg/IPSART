% PolC3 - 3x3 极化协方差矩阵
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef PolC3 < PolM3

    methods (Static)

        function C3 = fromS2(hh, hv, vh, vv)
            xx = (hv + vh) / sqrt(2);
            c11 = hh .* conj(hh);
            c12 = hh .* conj(xx);
            c13 = hh .* conj(vv);
            c22 = xx .* conj(xx);
            c23 = xx .* conj(vv);
            c33 = vv .* conj(vv);

            C3 = PolC3(c11, c22, c33, real(c12), real(c13), real(c23), ...
                imag(c12), imag(c13), imag(c23));
        end

    end

    methods

        function T3 = toT3(obj)
            C = zeros(3, 3, obj.Height, obj.Width, obj.Dtype);
            C(1,1,:,:) = obj.m11;
            C(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            C(1,3,:,:) = obj.m13_r + 1i * obj.m13_i;
            C(2,1,:,:) = conj(C(1,2,:,:));
            C(2,2,:,:) = obj.m22;
            C(2,3,:,:) = obj.m23_r + 1i * obj.m23_i;
            C(3,1,:,:) = conj(C(1,3,:,:));
            C(3,2,:,:) = conj(C(2,3,:,:));
            C(3,3,:,:) = obj.m33;

            M = cast([1 0 1; 1 0 -1; 0 sqrt(2) 0] / sqrt(2), obj.Dtype);
            T = pagemtimes(pagemtimes(M, C), M');

            T3 = PolT3(squeeze(T(1,1,:,:)), squeeze(T(2,2,:,:)), ...
                squeeze(T(3,3,:,:)), real(squeeze(T(1,2,:,:))), ...
                real(squeeze(T(1,3,:,:))), real(squeeze(T(2,3,:,:))), ...
                imag(squeeze(T(1,2,:,:))), imag(squeeze(T(1,3,:,:))), ...
                imag(squeeze(T(2,3,:,:))));
        end

    end

end