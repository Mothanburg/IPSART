% OrbitPolynomial - 轨道多项式对象
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
classdef OrbitPolynomial

    properties
        Coeff     % 轨道多项式系数
        Degree    % 轨道多项式阶数
        TimeRange % 轨道的时间范围
    end


    methods

        function obj = OrbitPolynomial(coefficent, degree, timeRange)
            arguments
                coefficent (:,6) double
                degree {mustBeInteger}
                timeRange (1,2) double
            end
            [rows,~] = size(coefficent);
            if rows - 1 ~= degree
                error("The degree must correspond to size of coefficient");
            end
            obj.Coeff = coefficent;
            obj.Degree = degree;
            obj.TimeRange = timeRange;
        end


        function T = TimeMatrix(obj, time, degree)
            arguments
                obj OrbitPolynomial
                time (:,1) double
                degree {mustBeInteger} = obj.Degree
            end

            t_norm = 4 * (time - obj.TimeRange(1)) / ...
                diff(obj.TimeRange) - 2;

            T = ones(length(t_norm), degree + 1);
            for n = 1:degree
                T(:,n + 1) = t_norm.^n;
            end
        end


        function state = GetSataliteState(obj, time, indices)
            arguments
                obj OrbitPolynomial
                time double
                indices (1,:) {mustBeInteger} = [1 2 3 4 5 6]
            end

            T = obj.TimeMatrix(time);

            state = T * obj.Coeff(:,indices);
        end


        function xyz = GetSatalitePosition(obj, relativeTime)
            xyz = obj.GetSataliteState(relativeTime, [1 2 3]);
        end


        function vxyz = GetSataliteVelocity(obj, relativeTime)
            vxyz = obj.GetSataliteState(relativeTime, [4 5 6]);
        end


    end


    methods (Static)

        % 根据轨道记录对象构造指定阶数（degree）的轨道多项式
        function obj = FromOrbitVector(orbitVector, degree)
            arguments
                orbitVector OrbitVector
                degree {mustBeInteger} = 5
            end

            % M = T * A
            % |X1 Y1 Z1 Vx1 Vy1 Vz1|   |1 t1 t1^2 ... t1^n|   |a0 b0 c0 d0 e0 f0|
            % |X2 Y2 Z2 Vx2 Vy2 Vz2|   |1 t2 t2^2 ... t2^n|   |a1 b1 c1 d1 e1 f1|
            % |X3 Y3 Z3 Vx3 Vy3 Vz3| = |1 t3 t3^2 ... t3^n| * |a2 b2 c2 d2 e2 f2|
            %          ...                     ...                    ...
            % |Xm Ym Zm Vxm Vym Vzm|   |1 tm tm^2 ... tm^n|   |an bn cn dn en fn|

            M = horzcat(orbitVector.X', orbitVector.Y', orbitVector.Z', ...
                orbitVector.Vx', orbitVector.Vy', orbitVector.Vz');

            t = orbitVector.Ticks';

            % rescale to [-2 2] for stability
            t = rescale(t, -2, 2);

            % generate T
            T = ones(orbitVector.Length, degree + 1);
            for n = 1:degree
                T(:,n + 1) = t.^n;
            end

            % calculate the matrix A
            obj = OrbitPolynomial(pinv(T) * M, degree, t([1 end]));
        end

    end

end
