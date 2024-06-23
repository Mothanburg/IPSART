classdef (Abstract) PolM2 < PolMat

    properties
        Height
        Width
    end

    properties (Access=protected)
        m11
        m12_r
        m12_i
        m22
    end

    properties (Dependent)
        Dtype
        SPAN
    end

    methods

        function obj = PolM2(m11, m12_r, m12_i, m22)
            [obj.Height,obj.Width] = size(m11);

            obj.m11 = zeros(obj.Height, obj.Width, class(m11));
            obj.m11(:,:) = m11;
            obj.m22 = zeros(obj.Height, obj.Width, class(m11));
            obj.m22(:,:) = m22;

            obj.m12_r = zeros(obj.Height, obj.Width, class(m11));
            obj.m12_r(:,:) = m12_r;
            obj.m12_i = zeros(obj.Height, obj.Width, class(m11));
            obj.m12_i(:,:) = m12_i;
        end

        function value = get.SPAN(obj)
            value = obj.m11 + obj.m22 + obj.m33;
        end

        function value = get.Dtype(obj)
            value = class(obj.m11);
        end

        function mat = getMatAt(obj, row, col)
            % assert(row >= 1 && row <= obj.Height && col >= 1 && col <= obj.Width);
            mat = zeros(2, obj.Dtype);
            mat(1,1) = obj.m11(row,col);
            mat(1,2) = obj.m12_r(row,col) + 1i * obj.m12_i(row,col);
            mat(2,1) = conj(mat(1,2));
            mat(2,2) = obj.m22(row,col);
        end

        function page = getPageAt(obj, x, y)
            % assert(x >= 1 && x <= 2 && y >= 1 && y <= 2);
            if x == y
                switch x
                    case 1
                        page = obj.m11;
                    case 2
                        page = obj.m22;
                end
            elseif x < y
                page = obj.m12_r + 1i * obj.m12_i;
            else
                page = conj(obj.getPageAt(y, x));
            end
        end

        function outObj = fmapPage(obj, func)
            outObj = obj;

            new_m11 = func(obj.m11);
            [outObj.Height,outObj.Width] = size(new_m11);

            outObj.m11 = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m11(:,:) = new_m11;

            new_m22 = func(obj.m22);
            outObj.m22 = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m22(:,:) = new_m22;

            new_m12 = func(obj.m12_r + 1i * obj.m12_i);
            outObj.m12_r = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m12_r(:,:) = real(new_m12);
            outObj.m12_i = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m12_i(:,:) = imag(new_m12);
        end

    end
end