classdef PolC3 < PolM3

    methods (Static)

        function C3 = fromS2(HH, HV, VH, VV)
            C3 = PolC3;

            [C3.Height,C3.Width] = size(HH);

            XX = (HV + VH) / sqrt(2);

            C3.m11 = HH .* conj(HH);
            C3.m22 = XX .* conj(XX);
            C3.m33 = VV .* conj(VV);

            C3.m12_r = real(HH) .* real(XX) + imag(HH) .* imag(XX);
            C3.m12_i = real(HH) .* -imag(XX) + real(XX) .* imag(HH);

            C3.m13_r = real(HH) .* real(VV) + imag(HH) .* imag(VV);
            C3.m13_i = real(HH) .* -imag(VV) + real(VV) .* imag(HH);

            C3.m23_r = real(XX) .* real(VV) + imag(XX) .* imag(VV);
            C3.m23_i = real(XX) .* -imag(VV) + real(VV) .* imag(XX);
        end
        
    end

    methods

        function T3 = toT3(obj)
            C(:,:,1,1) = obj.m11;
            C(:,:,1,2) = obj.m12_r + 1i * obj.m12_i;
            C(:,:,1,3) = obj.m13_r + 1i * obj.m13_i;
            C(:,:,2,1) = conj(C(:,:,1,2));
            C(:,:,2,2) = obj.m22;
            C(:,:,2,3) = obj.m23_r + 1i * obj.m23_i;
            C(:,:,3,1) = conj(C(:,:,1,3));
            C(:,:,3,2) = conj(C(:,:,2,3));
            C(:,:,3,3) = obj.m33;

            M = [1 0 1; 1 0 -1; 0 sqrt(2) 0] / sqrt(2);
            T = pagemtimes(pagemtimes(M, shiftdim(C, 2)), M');
            T = shiftdim(T, 2);

            T3 = PolT3.fromRaw(T(:,:,1,1), T(:,:,1,2), T(:,:,1,3), ...
                T(:,:,2,2), T(:,:,2,3), T(:,:,3,3));
        end

    end

end