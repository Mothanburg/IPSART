classdef SARDataReader

    properties (Abstract)
        FieldNames
    end

    properties (Dependent)
        NumOfLines    {mustBeInteger}
        NumOfPixels   {mustBeInteger}
        
        BandNumber
        BandInfo    (1,:) string

        Orbit             OrbitVector
        
        RangeSamplingRate double
        PRF               double

        CentralLongitude  double
        CentralLatitude   double
        SceneAltitude     double

        NearRange         double
        SceneStartTime    datetime

    end

    methods (Abstract)

        value = Query(obj, fieldName);
        data = ReadData(obj);

    end

    methods

        function value = get.BandNumber(obj)
            value = obj.Query("BandNumber");
        end
   
        function value = get.BandInfo(obj)
            value = obj.Query("BandInfo");
        end
    
    end

end