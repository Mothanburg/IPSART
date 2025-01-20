function clustered = internal__UnsupervisedWishartClustering_matlab( ...
    C, initClasses, numClasses, threshold, maxIter)

height = C.Height;
width = C.Width;

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
        c = C.Value.MatAt(i, j);
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


function vms = calc_centers(C, num_classes, classes)

vms = cell(1, num_classes);

C_full = zeros(C.Height, C.Width, C.Dim, C.Dim, C.Dtype);
for i = 1:C.Dim
    for j = 1:C.Dim
        C_full(:,:,i,j) = C.PageAt(i, j);
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