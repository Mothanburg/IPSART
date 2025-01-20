function [clustered, initClasses] = UnsupervisedWishartClustering(M, initClasses, options)

arguments
    M PolMat
    initClasses {mustBeInteger} = []
    options.THRESHOLD = 0.1
    options.MAX_ITER = inf
end

if isa(M, "PolC3") || isa(M, "PolC2")
    C = M;
elseif isa(M, "PolT3")
    C = M.toC3();
elseif isa(M, "PolT2")
    C = M.toC2();
else
    error("Input must can be converted to covariance matrix")
end

height = M.Height;
width = M.Width;

if isempty(initClasses)
    if isa(M, "PolC3")
        T3 = M.toT3();
    elseif isa(M, "PolT3")
        T3 = M;
    else
        error("You must specify an initial class if input is not 3x3 C or T")
    end

    num_classes = 8;
    [H,a,~] = CloudePottier(T3);
    initClasses = zeros(height, width, "int32");
    for i = 1:height
        for j = 1:width
            if H(i,j) <= 0.5
                if a(i,j) >= 47.5
                    initClasses(i,j) = 1;
                elseif a(i,j) >= 42.5
                    initClasses(i,j) = 2;
                else
                    initClasses(i,j) = 3;
                end
            elseif H(i,j) <= 0.9
                if a(i,j) >= 50
                    initClasses(i,j) = 4;
                elseif a(i,j) >= 40
                    initClasses(i,j) = 5;
                else
                    initClasses(i,j) = 6;
                end
            else
                if a(i,j) >= 55
                    initClasses(i,j) = 7;
                elseif a(i,j) >= 40
                    initClasses(i,j) = 8;
                else
                    % init_classes(i,j) = 9; % shouldn't happen
                    warning("The result of Cloude-Pottier has some problem")
                end
            end
        end
    end
else
    num_classes = numel(unique(initClasses(initClasses > 0)));
end

clustered = internal__UnsupervisedWishartClustering_matlab( ...
    C, initClasses, num_classes, options.THRESHOLD, options.MAX_ITER);

end
