classdef PolT3 < PolM3

    methods (Static)
        function T3 = fromS2(HH, HV, VH, VV)
            T3 = PolT3;

            [T3.Height,T3.Width] = size(HH);

            XX = (HV + VH) / sqrt(2);

            T3.m11 = (HH + VV) .* conj(HH + VV) / 2;
            T3.m22 = (HH - VV) .* conj(HH - VV) / 2;
            T3.m33 = 2 * XX .* conj(XX);

            t12 = (HH + VV) .* conj(HH - VV) / 2;
            T3.m12_r = zeros(T3.height, T3.width);
            T3.m12_r(:,:) = real(t12);
            T3.m12_i = zeros(T3.height, T3.width);
            T3.m12_i(:,:) = imag(t12);

            t13 = (HH + VV) .* conj(XX);
            T3.m13_r = zeros(T3.height, T3.width);
            T3.m13_r(:,:) = real(t12);
            T3.m13_i = zeros(T3.height, T3.width);
            T3.m13_i(:,:) = imag(t13);

            t23 = (HH - VV) .* conj(XX);
            T3.m23_r = zeros(T3.height, T3.width);
            T3.m23_r = real(t23);
            T3.m23_i = zeros(T3.height, T3.width);
            T3.m23_i = imag(t23);
        end
    end

    methods

        function C3 = toC3(obj)
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

            C3 = PolT3.fromRaw(C(:,:,1,1), C(:,:,1,2), C(:,:,1,3), ...
                C(:,:,2,2), C(:,:,2,3), C(:,:,3,3));
        end

    end

end