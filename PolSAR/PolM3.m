% PolM3 - 3x3 极化矩阵的抽象模型
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef (Abstract) PolM3 < PolMat

    properties
        Height
        Width
        Dim
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
        function obj = PolM3(m11, m22, m33, m12_r, m13_r, m23_r, m12_i, m13_i, m23_i)
            arguments
                m11 {mustBeReal}
                m22 {mustBeReal}
                m33 {mustBeReal}
                m12_r {mustBeReal}
                m13_r {mustBeReal}
                m23_r {mustBeReal}
                m12_i {mustBeReal}
                m13_i {mustBeReal}
                m23_i {mustBeReal}
            end
            obj.Dim = 3;
            [obj.Height,obj.Width] = size(m11);

            obj.m11 = m11;
            obj.m22 = m22;
            obj.m33 = m33;
            obj.m12_r = m12_r;
            obj.m13_r = m13_r;
            obj.m23_r = m23_r;
            obj.m12_i = m12_i;
            obj.m13_i = m13_i;
            obj.m23_i = m23_i;
        end

        function value = get.SPAN(obj)
            value = obj.m11 + obj.m22 + obj.m33;
        end

        function value = get.Dtype(obj)
            value = class(obj.m11);
        end

        function mat = MatAt(obj, row, col)
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

        function page = PageAt(obj, x, y)
            arguments
                obj
                x {mustBeInRange(x, 1, 3)}
                y {mustBeInRange(y, 1, 3)}
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
                page = conj(obj.PageAt(y, x));
            end
        end

        function outObj = MapPage(obj, func, varargin)
            cls = class(obj);

            new_m11 = func(obj.m11, varargin{:});
            new_m22 = func(obj.m22, varargin{:});
            new_m33 = func(obj.m33, varargin{:});
            new_m12_r = func(obj.m12_r, varargin{:});
            new_m12_i = func(obj.m12_i, varargin{:});
            new_m13_r = func(obj.m13_r, varargin{:});
            new_m13_i = func(obj.m13_i, varargin{:});
            new_m23_r = func(obj.m23_r, varargin{:});
            new_m23_i = func(obj.m23_i, varargin{:});

            outObj = feval(cls, new_m11, new_m22, new_m33, new_m12_r, ...
                new_m13_r, new_m23_r, new_m12_i, new_m13_i, new_m23_i);
        end

    end


end