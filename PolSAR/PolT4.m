classdef PolT4 < PolM4

    methods (Static)
        function T4 = fromS2(hh, hv, vh, vv)
            k1 = hh + vv;
            k2 = hh - vv;
            k3 = hv + vh;
            k4 = 1i * (hv - vh);

            t11 = k1 .* conj(k1);
            t12 = k1 .* conj(k2);
            t13 = k1 .* conj(k3);
            t14 = k1 .* conj(k4);
            t22 = k2 .* conj(k2);
            t23 = k2 .* conj(k3);
            t24 = k2 .* conj(k4);
            t33 = k3 .* conj(k3);
            t34 = k3 .* conj(k4);
            t44 = k4 .* conj(k4);

            T4 = PolT4(t11, t22, t33, t44, real(t12), real(t13), real(t14), ...
                real(t23), real(t24), real(t34), imag(t12), imag(t13), ...
                imag(t14), imag(t23), imag(t24), imag(t34));
        end
    end

    methods

        function C4 = toC4(obj)
            T = zeros(4, 4, obj.Height, obj.Width, obj.Dtype);
            T(1,1,:,:) = obj.m11;
            T(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            T(1,3,:,:) = obj.m13_r + 1i * obj.m13_i;
            T(1,4,:,:) = obj.m14_r + 1i * obj.m14_i;
            T(2,1,:,:) = conj(T(1,2,:,:));
            T(2,2,:,:) = obj.m22;
            T(2,3,:,:) = obj.m23_r + 1i * obj.m23_i;
            T(2,4,:,:) = obj.m24_r + 1i * obj.m24_i;
            T(3,1,:,:) = conj(T(1,3,:,:));
            T(3,2,:,:) = conj(T(2,3,:,:));
            T(3,3,:,:) = obj.m33;
            T(3,4,:,:) = obj.m34_r + 1i * obj.m34_i;
            T(4,1,:,:) = conj(T(1,4,:,:));
            T(4,2,:,:) = conj(T(2,4,:,:));
            T(4,3,:,:) = conj(T(3,4,:,:));
            T(4,4,:,:) = obj.m44;

            M = cast([1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]' / sqrt(2), obj.Dtype);
            C = pagemtimes(pagemtimes(M, T), M');
            clear T;

            C4 = PolC4(squeeze(C(1,1,:,:)), squeeze(C(2,2,:,:)), ...
                squeeze(C(3,3,:,:)), squeeze(C(4,4,:,:)),...
                real(squeeze(C(1,2,:,:))), real(squeeze(C(1,3,:,:))), ...
                real(squeeze(C(1,4,:,:))), real(squeeze(C(2,3,:,:))), ...
                real(squeeze(C(2,4,:,:))), real(squeeze(C(3,4,:,:))),...
                imag(squeeze(C(1,2,:,:))), imag(squeeze(C(1,3,:,:))), ...
                imag(squeeze(C(1,4,:,:))), imag(squeeze(C(2,3,:,:))), ...
                imag(squeeze(C(2,4,:,:))), imag(squeeze(C(3,4,:,:))));
        end

    end

end