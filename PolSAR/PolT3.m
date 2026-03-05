classdef PolT3 < PolM3

    methods (Static)
        function T3 = fromS2(hh, hv, vh, vv)
            xx = (hv + vh) / sqrt(2);
            t11 = (hh + vv) .* conj(hh + vv) / 2;
            t12 = (hh + vv) .* conj(hh - vv) / 2;
            t13 = (hh + vv) .* conj(xx);
            t22 = (hh - vv) .* conj(hh - vv) / 2;
            t23 = (hh - vv) .* conj(xx);
            t33 = 2 * xx .* conj(xx);

            T3 = PolT3(t11, t22, t33, real(t12), real(t13), real(t23), ...
                imag(t12), imag(t13), imag(t23));
        end
    end

    methods

        function C3 = toC3(obj)
            T = zeros(3, 3, obj.Height, obj.Width, obj.Dtype);
            T(1,1,:,:) = obj.m11;
            T(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            T(1,3,:,:) = obj.m13_r + 1i * obj.m13_i;
            T(2,1,:,:) = conj(T(1,2,:,:));
            T(2,2,:,:) = obj.m22;
            T(2,3,:,:) = obj.m23_r + 1i * obj.m23_i;
            T(3,1,:,:) = conj(T(1,3,:,:));
            T(3,2,:,:) = conj(T(2,3,:,:));
            T(3,3,:,:) = obj.m33;

            M = cast([1 0 1; 1 0 -1; 0 sqrt(2) 0]' / sqrt(2), obj.Dtype);
            C = pagemtimes(pagemtimes(M, T), M');

            C3 = PolC3(squeeze(C(1,1,:,:)), squeeze(C(2,2,:,:)), ...
                squeeze(C(3,3,:,:)), real(squeeze(C(1,2,:,:))), ...
                real(squeeze(C(1,3,:,:))), real(squeeze(C(2,3,:,:))), ...
                imag(squeeze(C(1,2,:,:))), imag(squeeze(C(1,3,:,:))), ...
                imag(squeeze(C(2,3,:,:))));
        end

    end

end