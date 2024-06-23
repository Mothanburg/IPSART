classdef (Abstract) PolMat

    properties (Abstract)
        Height
        Width
    end

    properties (Abstract, Dependent)
        Dtype
        SPAN
    end

    methods (Abstract)
        mat = getMatAt(obj, row, col)
        page = getPageAt(obj, x, y)
        outObj = fmapPage(obj, func)
    end

end