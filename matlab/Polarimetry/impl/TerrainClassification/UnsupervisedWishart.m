% 无监督Wishart分类算法
function classes = UnsupervisedWishart(C, numClasses, initClasses, options)

arguments
    C
    numClasses
    initClasses
    options.THRESHOLD = 0.1
    options.MAX_ITER = inf
end

[height,width,~,~] = size(C);

transfer_rate = 1;
iter_input = initClasses;
i = 1;
len = height*width;
while transfer_rate > options.THRESHOLD && i <= options.MAX_ITER
    vms = calc_centers(numClasses, C, iter_input);

    iter_output = zeros(height, width);
    for j=1:len
        [m,n] = ind2sub([height, width], j);
        c = squeeze(C(m,n,:,:));
        distances = calc_distances(numClasses, c, vms);
        [~,new_class] = min(distances);
        iter_output(j) = new_class;
    end
    transfer_rate = sum(iter_output~=iter_input, 'all')/(height*width);

    iter_input = iter_output;
    i = i+1;
end
classes = iter_output;

end


% 计算类别中心
function vms = calc_centers(num_classes, C, classes)

[~,~,dim,~] = size(C);
vms = cell(1, num_classes);
for m = 1:num_classes
    vm = zeros(dim, dim);
    mask = classes==m;
    count = sum(mask, "all");
    if count ~= 0
        C_masked = C.*mask;
        vm = squeeze(sum(C_masked, [1 2])/count);
    end
    vms{m} = vm;
end

end


% 计算与各类别的距离
function ds = calc_distances(num_classes, c, vms)

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

