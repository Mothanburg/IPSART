% OrbitVector - 轨道记录对象
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef OrbitVector

    properties
        X (1,:) double
        Y (1,:) double
        Z (1,:) double
        Vx (1,:) double
        Vy (1,:) double
        Vz (1,:) double
        Ticks (1,:) double
        TimeStamp (1,:) datetime
        EcefEllipsoid referenceEllipsoid
    end

    properties (Dependent)
        Length
    end

    methods


        function obj = OrbitVector(x, y, z, vx, vy, vz, timestamp, reference)
            len = length(x);
            assert(length(y) == len && ...
                   length(z) == len && ...
                   length(vx) == len && ...
                   length(vy) == len && ...
                   length(vz) == len && ...
                   length(timestamp) == len, ...
                "The length of input paramters must be same");

            % sort the records in ascending order of timestamp
            [obj.TimeStamp,indices] = sort(timestamp);
            obj.X = x(indices);
            obj.Y = y(indices);
            obj.Z = z(indices);
            obj.Vx = vx(indices);
            obj.Vy = vy(indices);
            obj.Vz = vz(indices);

            % calculate relative time
            time_diff = obj.TimeStamp - obj.TimeStamp(1);
            obj.Ticks = time_diff.Second;
            obj.EcefEllipsoid = reference;
        end


        function value = get.Length(obj)
            value = length(obj.X);
        end


        function tick = TickOf(obj, timestamp)
            arguments
                obj OrbitVector
                timestamp datetime
            end

            time_diff = timestamp - obj.TimeStamp(1);
            tick = time_diff.Second;
        end


    end

end