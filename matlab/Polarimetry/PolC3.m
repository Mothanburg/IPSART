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

            C3 = PolC3(c11, real(c12), imag(c12), real(c13), imag(c13), c22, ...
                real(c23), imag(c23), c33);
        end

    end

    methods

        function T3 = toT3(obj)
            C = zeros(obj.Height, obj.Width, 3, 3, obj.Dtype);
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

            T3 = PolT3(T(:,:,1,1), real(T(:,:,1,2)), imag(T(:,:,1,2)), ...
                real(T(:,:,1,3)), imag(T(:,:,1,3)), T(:,:,2,2), real(T(:,:,2,3)), ...
                imag(T(:,:,2,3)), T(:,:,3,3));
        end

    end

end