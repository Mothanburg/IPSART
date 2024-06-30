classdef PolT2 < PolM2

    methods (Static)

        function T2 = fromHHVV(hh, vv)
            t11 = (hh + vv) .* conj(hh + vv);
            t12 = (hh + vv) .* conj(hh - vv);
            t22 = (hh - vv) .* conj(hh - vv);

            T2 = PolT2(t11, t22, real(t12), imag(t12));
        end

    end

    methods

        function C2 = toC2(obj)
            T = zeros(obj.Height, obj.Width, 2, 2, obj.Dtype);
            T(:,:,1,1) = obj.m11;
            T(:,:,1,2) = obj.m12_r + 1i * obj.m12_i;
            T(:,:,2,1) = conj(T(:,:,1,2));
            T(:,:,2,2) = obj.m22;
            M = [1, 1; 1, -1]';
            C = pagemtimes(pagemtimes(M, shiftdim(T, 2)), M');
            C = shiftdim(C, 2);

            C2 = PolC2(C(:,:,1,1), C(:,:,2,2), real(C(:,:,1,2)), ...
                imag(C(:,:,1,2)), "HHVV");
        end

    end

end