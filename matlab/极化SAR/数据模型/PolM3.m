classdef PolM3 < PolMat

    properties (Access = protected)
        m11 double
        m22 double
        m33 double
        m12_r double
        m13_r double
        m23_r double
        m12_i double
        m13_i double
        m23_i double
    end

    properties (Dependent)
        SPAN double
    end

    methods (Static)
        function M3 = fromRaw(m11, m12, m13, m22, m23, m33)
            M3 = PolM3;
            [M3.Height,M3.Width] = size(m11);
            M3.m11 = m11;
            M3.m22 = m22;
            M3.m33 = m33;

            M3.m12_r = zeros(M3.Height,M3.Width);
            M3.m12_r(:,:) = real(m12);
            M3.m13_r = zeros(M3.Height,M3.Width);
            M3.m13_r(:,:) = real(m13);
            M3.m23_r = zeros(M3.Height,M3.Width);
            M3.m23_r(:,:) = real(m23);

            M3.m12_i = zeros(M3.Height,M3.Width);
            M3.m12_i(:,:) = imag(m12);
            M3.m13_i = zeros(M3.Height,M3.Width);
            M3.m13_i(:,:) = imag(m13);
            M3.m23_i = zeros(M3.Height,M3.Width);
            M3.m23_i(:,:) = imag(m23);
        end
    end


    methods
        function value = get.SPAN(obj)
            value = obj.m11 + obj.m22 + obj.m33;
        end

        function mat = getMatAt(obj, row, col)
            mat = zeros(3);
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
            arguments
                obj PolM3
                x int {mustBeLessThanOrEqual(x, 3), mustBeGreaterThanOrEqual(x, 1)}
                y int {mustBeLessThanOrEqual(y, 3), mustBeGreaterThanOrEqual(y, 1)}
            end

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
            outObj = PolM3;
            outObj.Height = obj.Height;
            outObj.Width = obj.Width;

            outObj.m11 = zeros(obj.Height, obj.Width);
            outObj.m11(:,:) = func(obj.m11);

            outObj.m22 = zeros(obj.Height, obj.Width);
            outObj.m22(:,:) = func(obj.m22);

            outObj.M33 = zeros(obj.Height, obj.Width);
            outObj.M33(:,:) = func(obj.m33);

            outObj.m12_r = zeros(obj.Height, obj.Width);
            outObj.m12_i = zeros(obj.Height, obj.Width);
            m12 = func(obj.m12_r + 1i * obj.m12_i);
            outObj.m12_r(:,:) = real(m12);
            outObj.m12_i(:,:) = imag(m12);

            outObj.m13_r = zeros(obj.Height, obj.Width);
            outObj.m13_i = zeros(obj.Height, obj.Width);
            m13 = func(obj.m13_r + 1i * obj.m13_i);
            outObj.m13_r(:,:) = real(m13);
            outObj.m13_i(:,:) = imag(m13);

            outObj.m23_r = zeros(obj.Height, obj.Width);
            outObj.m23_i = zeros(obj.Height, obj.Width);
            m23 = func(obj.m23_r + 1i * obj.m23_i);
            outObj.m23_r(:,:) = real(m23);
            outObj.m23_i(:,:) = imag(m23);
        end

    end


end