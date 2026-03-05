classdef PolC2 < PolM2

    properties
        PolType string
    end

    methods (Static)

        function C2 = fromChannel(Ch1, Ch2, polType)
            arguments
                Ch1
                Ch2
                polType string = "Unknown"
            end
            c11 = Ch1 .* conj(Ch1);
            c12 = Ch1 .* conj(Ch2);
            c22 = Ch2 .* conj(Ch2);
            C2 = PolC2(c11, c22, real(c12), imag(c12), polType);
        end

    end

    methods

        function obj = PolC2(m11, m22, m12_r, m12_i, polType)
            obj = obj@PolM2(m11, m22, m12_r, m12_i);
            obj.PolType = upper(polType);
        end

        function T2 = toT2(obj)
            if obj.PolType ~= "HHVV"
                error("Only ""HHVV"" is compitable with the T matrix.")
            end
            C = zeros(2, 2, obj.Height, obj.Width, obj.Dtype);
            C(1,1,:,:) = obj.m11;
            C(1,2,:,:) = obj.m12_r + 1i * obj.m12_i;
            C(2,1,:,:) = conj(C(1,2,:,:));
            C(2,2,:,:) = obj.m22;

            M = cast([1, 1; 1, -1], obj.Dtype);
            T = pagemtimes(pagemtimes(M, C), M');

            T2 = PolT2(squeeze(T(1,1,:,:)), squeeze(T(2,2,:,:)), ...
                real(squeeze(T(1,2,:,:))), imag(squeeze(T(1,2,:,:))));
        end

    end

end