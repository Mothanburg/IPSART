% PolC4 - 4x4 极化协方差矩阵
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef PolC4 < PolM4

    methods (Static)

        function C4 = fromS2(hh, hv, vh, vv)
            c11 = hh .* conj(hh);
            c12 = hh .* conj(hv);
            c13 = hh .* conj(vh);
            c14 = hh .* conj(vv);
            c22 = hv .* conj(hv);
            c23 = hv .* conj(vh);
            c24 = hv .* conj(vv);
            c33 = vh .* conj(vh);
            c34 = vh .* conj(vv);
            c44 = vv .* conj(vv);

            C4 = PolC4(c11, c22, c33, c44, real(c12), real(c13), real(c14), ...
                real(c23), real(c24), real(c34), imag(c12), imag(c13), ...
                imag(c14), imag(c23), imag(c24), imag(c34));
        end

    end

    methods

        function T4 = toT4(obj)
            C = zeros(4, 4, obj.Height, obj.Width, obj.Dtype);
            C(1,1,:,:) = obj.m11;
            C(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            C(1,3,:,:) = obj.m13_r + 1i * obj.m13_i;
            C(1,4,:,:) = obj.m14_r + 1i * obj.m14_i;
            C(2,1,:,:) = conj(C(1,2,:,:));
            C(2,2,:,:) = obj.m22;
            C(2,3,:,:) = obj.m23_r + 1i * obj.m23_i;
            C(2,4,:,:) = obj.m24_r + 1i * obj.m24_i;
            C(3,1,:,:) = conj(C(1,3,:,:));
            C(3,2,:,:) = conj(C(2,3,:,:));
            C(3,3,:,:) = obj.m33;
            C(3,4,:,:) = obj.m34_r + 1i * obj.m34_i;
            C(4,1,:,:) = conj(C(1,4,:,:));
            C(4,2,:,:) = conj(C(2,4,:,:));
            C(4,3,:,:) = conj(C(3,4,:,:));
            C(4,4,:,:) = obj.m44;

            M = cast([1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0] / sqrt(2), obj.Dtype);
            T = pagemtimes(pagemtimes(M, C), M');
            clear C;

            T4 = PolT4(squeeze(T(1,1,:,:)), squeeze(T(2,2,:,:)), ...
                squeeze(T(3,3,:,:)), squeeze(T(4,4,:,:)),...
                real(squeeze(T(1,2,:,:))), real(squeeze(T(1,3,:,:))), ...
                real(squeeze(T(1,4,:,:))), real(squeeze(T(2,3,:,:))), ...
                real(squeeze(T(2,4,:,:))), real(squeeze(T(3,4,:,:))),...
                imag(squeeze(T(1,2,:,:))), imag(squeeze(T(1,3,:,:))), ...
                imag(squeeze(T(1,4,:,:))), imag(squeeze(T(2,3,:,:))), ...
                imag(squeeze(T(2,4,:,:))), imag(squeeze(T(3,4,:,:))));
        end

    end

end