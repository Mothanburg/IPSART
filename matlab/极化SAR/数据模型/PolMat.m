classdef (Abstract) PolMat

    properties (Abstract)
        Height int
        Width int
    end

    methods
        span = getSPAN(obj)
        mat = getMatAt(obj, row, col)
        page = getPageAt(obj, x, y)
        outObj = fmapPage(obj, func)
    end

end