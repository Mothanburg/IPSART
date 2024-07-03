function [clustered, initClasses] = UnsupervisedWishartClustering(C, initClasses, options)

arguments
    C PolMat
    initClasses {mustBeInteger} = []
    options.THRESHOLD = 0.1
    options.MAX_ITER = inf
end

if ~isa(C, "PolC3") && ~isa(C, "PolC2")
    error("Only covariance matrix C can apply the Wishart classifier.")
end

height = C.Height;
width = C.Width;

if isempty(initClasses)
    if isa(C, "PolC3")
        num_classes = 8;
        [H,a,~] = CloudePottier(C.toT3());
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
                        % init_classes(i,j) = 9; % 不存在
                        warning("The result of Cloude-Pottier has some problem.")
                    end
                end
            end
        end
    else
        error("You must specify an initial class if input is not 3x3 covariance matrix C.")
    end
else
    num_classes = numel(unique(initClasses(initClasses > 0)));
end

clustered = UWC_matlab(C, initClasses, num_classes, options.THRESHOLD, options.MAX_ITER);

end


function clustered = UWC_matlab(C, initClasses, numClasses, threshold, maxIter)

transfer_rate = 1;
iter_input = initClasses;

iters = 1;
len = height * width;
C = parallel.pool.Constant(C);
while transfer_rate > threshold && iters <= maxIter
    vms = calc_centers(C.Value, numClasses, iter_input);

    iter_output = zeros(height, width, "int32");
    parfor idx = 1:len
        [i,j] = ind2sub([height, width], idx);
        c = C.Value.getMatAt(i, j);
        distances = calc_distances(c, vms);
        [~,new_class] = min(distances);
        iter_output(idx) = new_class;
    end
    transfer_rate = sum(iter_output~=iter_input, 'all') / len;

    iter_input = iter_output;
    iters = iters + 1;
end
clustered = iter_output;

end

% 计算类别中心
function vms = calc_centers(C, num_classes, classes)

vms = cell(1, num_classes);

C_full = zeros(C.Height, C.Width, C.Dim, C.Dim, C.Dtype);
for i = 1:C.Dim
    for j = 1:C.Dim
        C_full(:,:,i,j) = C.getPageAt(i, j);
    end
end

for m = 1:num_classes
    vm = zeros(C.Dim, C.Dim, C.Dtype);
    mask = classes==m;
    count = sum(mask, "all");
    if count ~= 0
        C_masked = C_full .* mask;
        vm = squeeze(sum(C_masked, [1 2])/count);
    end
    vms{m} = vm;
end

end


% 计算与各类别的距离
function ds = calc_distances(c, vms)

num_classes = length(vms);
ds = zeros(1, num_classes);

for m = 1:num_classes
    vm = vms{m};
    if all(vm==0, "all")
        ds(m) = inf;
    else
        ds(m) = log(real(det(vm)))+real(trace(vm\c));
    end
end

end

