classdef (Abstract) PolM2 < PolMat

    properties
        Height
        Width
        Dim
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

        function obj = PolM2(m11, m22, m12_r, m12_i)
            obj.Dim = 2;
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
            value = obj.m11 + obj.m22;
        end

        function value = get.Dtype(obj)
            value = class(obj.m11);
        end

        function mat = MatAt(obj, row, col)
            assert(row >= 1 && row <= obj.Height && col >= 1 && col <= obj.Width);
            mat = zeros(2, obj.Dtype);
            mat(1,1) = obj.m11(row,col);
            mat(1,2) = obj.m12_r(row,col) + 1i * obj.m12_i(row,col);
            mat(2,1) = conj(mat(1,2));
            mat(2,2) = obj.m22(row,col);
        end

        function page = PageAt(obj, x, y)
            assert(x >= 1 && x <= 2 && y >= 1 && y <= 2);
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
                page = conj(obj.PageAt(y, x));
            end
        end

        function outObj = MapPage(obj, func, varargin)
            cls = class(obj);

            new_m11 = func(obj.m11, varargin{:});
            new_m22 = func(obj.m22, varargin{:});
            new_m12_r = func(obj.m12_r, varargin{:});
            new_m12_i = func(obj.m12_i, varargin{:});

            if cls == "PolC2"
                outObj = feval(cls, new_m11, new_m22, new_m12_r, new_m12_i, obj.PolType);
            else
                outObj = feval(cls, new_m11, new_m22, new_m12_r, new_m12_i);
            end
        end

    end
end