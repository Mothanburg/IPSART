classdef PolC2 < PolM2

    properties
        PolType string
    end

    methods (Static)

        function C2 = fromPolChan(polCh1, polCh2, polType)
            c11 = polCh1 .* conj(polCh1);
            c12 = polCh1 .* conj(polCh2);
            c22 = polCh2 .* conj(polCh2);
            C2 = PolC2(c11, real(c12), imag(c12), c22, polType);
        end

    end

    methods

        function obj = PolC2(m11, m12, m22, polType)
            obj = obj@PolM2(m11, m12, m22);
            obj.PolType = polType;
        end

        function T2 = toT2(obj)
            if obj.PolType ~= "HHVV"
                error("Only ""HHVV"" is compitable with the T matrix.")
            end
            C = zeros(obj.Height, obj.Width, 2, 2, obj.Dtype);
            C(:,:,1,1) = obj.m11;
            C(:,:,1,2) = obj.m12_r + 1i * obj.m12_i;
            C(:,:,2,1) = conj(C(:,:,1,2));
            C(:,:,2,2) = obj.m22;
            M = [1, 1; 1, -1];
            T = pagemtimes(pagemtimes(M, shiftdim(C, 2)), M');
            T = shiftdim(T, 2);

            T2 = PolT2(T(:,:,1,1), real(T(:,:,1,2)), imag(T(:,:,1,2), ...
                T(:,:,2,2));
        end

    end

end