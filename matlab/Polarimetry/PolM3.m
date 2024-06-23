classdef (Abstract) PolM3 < PolMat

    properties
        Height
        Width
    end

    properties (Access=protected)
        m11
        m22
        m33
        m12_r
        m13_r
        m23_r
        m12_i
        m13_i
        m23_i
    end

    properties (Dependent)
        Dtype
        SPAN
    end

    methods
        function obj = PolM3(m11, m12_r, m12_i, m13_r, m13_i, m22, m23_r, m23_i, m33)
            [obj.Height,obj.Width] = size(m11);

            obj.m11 = zeros(obj.Height, obj.Width, class(m11));
            obj.m11(:,:) = real(m11);
            obj.m22 = zeros(obj.Height, obj.Width, class(m11));
            obj.m22(:,:) = real(m22);
            obj.m33 = zeros(obj.Height, obj.Width, class(m11));
            obj.m33(:,:) = real(m33);

            obj.m12_r = zeros(obj.Height, obj.Width, class(m11));
            obj.m12_r(:,:) = m12_r;
            obj.m13_r = zeros(obj.Height, obj.Width, class(m11));
            obj.m13_r(:,:) = m13_r;
            obj.m23_r = zeros(obj.Height, obj.Width, class(m11));
            obj.m23_r(:,:) = m23_r;

            obj.m12_i = zeros(obj.Height, obj.Width, class(m11));
            obj.m12_i(:,:) = m12_i;
            obj.m13_i = zeros(obj.Height, obj.Width, class(m11));
            obj.m13_i(:,:) = m13_i;
            obj.m23_i = zeros(obj.Height, obj.Width, class(m11));
            obj.m23_i(:,:) = m23_i;
        end

        function value = get.SPAN(obj)
            value = obj.m11 + obj.m22 + obj.m33;
        end

        function value = get.Dtype(obj)
            % assert(all(class(obj.m11) == class(obj.m23_r)) * ...
            %     all(class(obj.m22) == class(obj.m13_i)));
            value = class(obj.m11);
        end

        function mat = getMatAt(obj, row, col)
            % assert(row >= 1 && row <= obj.Height && col >= 1 && col <= obj.Width);
            mat = zeros(3, obj.Dtype);
            mat(1,1) = obj.m11(row,col);
            mat(1,2) = obj.m12_r(row,col) + 1i * obj.m12_i(row,col);
            mat(1,3) = obj.m13_r(row,col) + 1i * obj.m13_i(row,col);
            mat(2,1) = conj(mat(1,2));
            mat(2,2) = obj.m22(row,col);
            mat(2,3) = obj.m23_r(row,col) + 1i * obj.m23_i(row,col);
            mat(3,1) = conj(mat(1,3));
            mat(3,2) = conj(mat(2,3));
            mat(3,3) = obj.m33(row,col);
        end

        function page = getPageAt(obj, x, y)
            % assert(x >= 1 && x <= 3 && y >= 1 && y <= 3);
            if x == y
                switch x
                    case 1
                        page = obj.m11;
                    case 2
                        page = obj.m22;
                    case 3
                        page = obj.m33;
                end
            elseif x < y
                if x == 1
                    if y == 2
                        page = obj.m12_r + 1i * obj.m12_i;
                    else % y == 3
                        page = obj.m13_r + 1i * obj.m13_i;
                    end
                else % x == 2, y == 3
                    page = obj.m23_r + 1i * obj.m23_i;
                end
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

            new_m33 = func(obj.m33);
            outObj.m33 = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m33(:,:) = new_m33;

            new_m12 = func(obj.m12_r + 1i * obj.m12_i);
            outObj.m12_r = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m12_r(:,:) = real(new_m12);
            outObj.m12_i = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m12_i(:,:) = imag(new_m12);

            new_m13 = func(obj.m13_r + 1i * obj.m13_i);
            outObj.m13_r = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m13_r(:,:) = real(new_m13);
            outObj.m13_i = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m13_i(:,:) = imag(new_m13);

            new_m23 = func(obj.m23_r + 1i * obj.m23_i);
            outObj.m23_r = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m23_r(:,:) = real(new_m23);
            outObj.m23_i = zeros(outObj.Height, outObj.Width, obj.Dtype);
            outObj.m23_i(:,:) = imag(new_m23);
        end

    end


end