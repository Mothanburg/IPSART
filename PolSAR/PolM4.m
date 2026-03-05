classdef (Abstract) PolM4 < PolMat

    properties
        Height
        Width
        Dim
        m11
        m22
        m33
        m44
        m12_r
        m13_r
        m14_r
        m23_r
        m24_r
        m34_r
        m12_i
        m13_i
        m14_i
        m23_i
        m24_i
        m34_i
    end

    properties (Dependent)
        Dtype
        SPAN
    end

    methods
        function obj = PolM3(m11, m22, m33, m44, m12_r, m13_r, m14_r, m23_r, ...
                m24_r, m34_r, m12_i, m13_i, m14_i, m23_i, m24_i, m34_i)
            arguments
                m11 {mustBeReal}
                m22 {mustBeReal}
                m33 {mustBeReal}
                m44 {mustBeReal}
                m12_r {mustBeReal}
                m13_r {mustBeReal}
                m14_r {mustBeReal}
                m23_r {mustBeReal}
                m24_r {mustBeReal}
                m34_r {mustBeReal}
                m12_i {mustBeReal}
                m13_i {mustBeReal}
                m14_i {mustBeReal}
                m23_i {mustBeReal}
                m24_i {mustBeReal}
                m34_i {mustBeReal}
            end
            obj.Dim = 4;
            [obj.Height,obj.Width] = size(m11);

            obj.m11 = m11;
            obj.m22 = m22;
            obj.m33 = m33;
            obj.m44 = m44;
            obj.m12_r = m12_r;
            obj.m13_r = m13_r;
            obj.m14_r = m14_r;
            obj.m23_r = m23_r;
            obj.m24_r = m24_r;
            obj.m34_r = m34_r;
            obj.m12_i = m12_i;
            obj.m13_i = m13_i;
            obj.m14_i = m14_i;
            obj.m23_i = m23_i;
            obj.m24_i = m24_i;
            obj.m34_i = m34_i;
        end

        function value = get.SPAN(obj)
            value = obj.m11 + obj.m22 + obj.m33 + obj.m44;
        end

        function value = get.Dtype(obj)
            value = class(obj.m11);
        end

        function mat = MatAt(obj, row, col)
            mat = zeros(4, obj.Dtype);
            mat(1,1) = obj.m11(row,col);
            mat(1,2) = obj.m12_r(row,col) + 1i * obj.m12_i(row,col);
            mat(1,3) = obj.m13_r(row,col) + 1i * obj.m13_i(row,col);
            mat(1,4) = obj.m14_r(row,col) + 1i * obj.m14_i(row,col);
            mat(2,1) = conj(mat(1,2));
            mat(2,2) = obj.m22(row,col);
            mat(2,3) = obj.m23_r(row,col) + 1i * obj.m23_i(row,col);
            mat(2,4) = obj.m24_r(row,col) + 1i * obj.m24_i(row,col);
            mat(3,1) = conj(mat(1,3));
            mat(3,2) = conj(mat(2,3));
            mat(3,3) = obj.m33(row,col);
            mat(3,4) = obj.m34_r(row,col) + 1i * obj.m34_i(row,col);
            mat(4,1) = conj(mat(1,4));
            mat(4,2) = conj(mat(2,4));
            mat(4,3) = conj(mat(3,4));
            mat(4,4) = obj.m44(row,col);
        end

        function page = PageAt(obj, x, y)
            arguments
                obj
                x {mustBeInRange(x, 1, 4)}
                y {mustBeInRange(y, 1, 4)}
            end

            if x == y
                switch x
                    case 1
                        page = obj.m11;
                    case 2
                        page = obj.m22;
                    case 3
                        page = obj.m33;
                    case 4
                        page = obj.m44;
                end
            elseif x < y
                if x == 1
                    switch y
                        case 2
                            page = obj.m12_r + 1i * obj.m12_i;
                        case 3
                            page = obj.m13_r + 1i * obj.m13_i;
                        case 4
                            page = obj.m14_r + 1i * obj.m14_i;
                    end
                elseif x == 2
                    if y == 3
                        page = obj.m23_r + 1i * obj.m23_i;
                    else % (x,y) = (2,4)
                        page = obj.m24_r + 1i * obj.m24_i;
                    end
                else % (x,y) = (3,4)
                    page = obj.m34_r + 1i * obj.m34_i;
                end
            else
                page = conj(obj.PageAt(y, x));
            end

        end

        function outObj = MapPage(obj, func, varargin)
            cls = class(obj);

            new_m11 = func(obj.m11, varargin{:});
            new_m22 = func(obj.m22, varargin{:});
            new_m33 = func(obj.m33, varargin{:});
            new_m44 = func(obj.m44, varargin{:});
            new_m12_r = func(obj.m12_r, varargin{:});
            new_m12_i = func(obj.m12_i, varargin{:});
            new_m13_r = func(obj.m13_r, varargin{:});
            new_m13_i = func(obj.m13_i, varargin{:});
            new_m14_r = func(obj.m14_r, varargin{:});
            new_m14_i = func(obj.m14_i, varargin{:});
            new_m23_r = func(obj.m23_r, varargin{:});
            new_m23_i = func(obj.m23_i, varargin{:});
            new_m24_r = func(obj.m24_r, varargin{:});
            new_m24_i = func(obj.m24_i, varargin{:});
            new_m34_r = func(obj.m34_r, varargin{:});
            new_m34_i = func(obj.m34_i, varargin{:});

            outObj = feval(cls, new_m11, new_m22, new_m33, new_m44, ...
                new_m12_r, new_m13_r, new_m14_r, new_m23_r, new_m24_r, ...
                new_m34_r, new_m12_i, new_m13_i, new_m14_i, new_m23_i, ...
                new_m24_i, new_m34_i);
        end

    end


end