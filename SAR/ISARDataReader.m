classdef ISARDataReader

    properties (Abstract)
        FieldNames
    end

    methods (Abstract)

        value = Query(obj, fieldName);
        data = ReadData(obj);

    end

end