classdef PolT3 < PolM3

    methods (Static)
        function T3 = fromS2(hh, hv, vh, vv)
            arguments (Input)
                hh (:,:)
                hv (:,:)
                vh (:,:)
                vv (:,:)
            end

            arguments (Output)
                T3 PolT3
            end

            xx = (hv + vh) / sqrt(2);
            t11 = (hh + vv) .* conj(hh + vv) / 2;
            t12 = (hh + vv) .* conj(hh - vv) / 2;
            t13 = (hh + vv) .* conj(xx);
            t22 = (hh - vv) .* conj(hh - vv) / 2;
            t23 = (hh - vv) .* conj(xx);
            t33 = 2 * xx .* conj(xx);

            T3 = PolT3(t11, real(t12), imag(t12), real(t13), imag(t13), t22, ...
                real(t23), imag(t23), t33);
        end
    end

    methods

        function C3 = toC3(obj)
            T = zeros(obj.Height, obj.Width, 3, 3, obj.Dtype);
            T(:,:,1,1) = obj.m11;
            T(:,:,1,2) = obj.m12_r + 1i * obj.m12_i;
            T(:,:,1,3) = obj.m13_r + 1i * obj.m13_i;
            T(:,:,2,1) = conj(T(:,:,1,2));
            T(:,:,2,2) = obj.m22;
            T(:,:,2,3) = obj.m23_r + 1i * obj.m23_i;
            T(:,:,3,1) = conj(T(:,:,1,3));
            T(:,:,3,2) = conj(T(:,:,2,3));
            T(:,:,3,3) = obj.m33;

            M = [1 0 1; 1 0 -1; 0 sqrt(2) 0]' / sqrt(2);
            C = pagemtimes(pagemtimes(M, shiftdim(T, 2)), M');
            C = shiftdim(C, 2);

            C3 = PolC3(C(:,:,1,1), real(C(:,:,1,2)), imag(C(:,:,1,2)), ...
                real(C(:,:,1,3)), imag(C(:,:,1,3)), C(:,:,2,2), real(C(:,:,2,3)), ...
                imag(C(:,:,2,3)), C(:,:,3,3));
        end

    end

end