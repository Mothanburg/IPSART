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
            T = zeros(2, 2, obj.Height, obj.Width, obj.Dtype);
            T(1,1,:,:) = obj.m11;
            T(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            T(2,1,:,:) = conj(T(1,2,:,:));
            T(2,2,:,:) = obj.m22;
            
            M = cast([1, 1; 1, -1]', obj.Dtype);
            C = pagemtimes(pagemtimes(M, T), M');

            C2 = PolC2(squeeze(C(1,1,:,:)), squeeze(C(2,2,:,:)), ...
                real(squeeze(C(1,2,:,:))), imag(squeeze(C(1,2,:,:))), "HHVV");
        end

    end

end