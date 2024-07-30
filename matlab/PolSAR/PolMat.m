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
        outObj = MapPage(obj, func)
    end

end