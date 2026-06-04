% PolMat - 极化矩阵的抽象模型
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef (Abstract) PolMat

    properties (Abstract)
        Height
        Width
        Dim
    end

    properties (Abstract, Dependent)
        Dtype
        SPAN
    end

    methods (Abstract)
        mat = MatAt(obj, row, col)
        page = PageAt(obj, x, y)
        outObj = MapPage(obj, func) % the func must be a linear operator
    end

end