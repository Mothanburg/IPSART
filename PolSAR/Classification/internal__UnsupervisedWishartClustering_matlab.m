% internal__UnsupervisedWishartClustering_matlab - Internal MATLAB implementation of Wishart clustering
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function clustered = internal__UnsupervisedWishartClustering_matlab( ...
    C, initClasses, numClasses, threshold, maxIter)

height = C.Height;
width = C.Width;

transfer_rate = 1;
iter_input = initClasses;

iters = 1;
len = height * width;
while transfer_rate > threshold && iters <= maxIter
    vms = zeros(C.Dim, C.Dim, numClasses);
    for i = 1:C.Dim
        for j = i:C.Dim
            page = C.PageAt(i,j);
            for cls = 1:numClasses
                vms(i,j,cls) = mean(page(iter_input == cls));
                vms(j,i,cls) = conj(vms(i,j,cls));
            end
        end
    end

    iter_output = zeros(height, width);
    parfor idx = 1:len
        [i,j] = ind2sub([height, width], idx);
        c = C.MatAt(i, j);
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


function ds = calc_distances(c, vms)
num_classes = size(vms, 3);
ds = zeros(1, num_classes);

for m = 1:num_classes
    vm = squeeze(vms(:,:,m));
    if any(isnan(vm(:)))
        ds(m) = inf;
    else
        ds(m) = log(real(det(vm)))+real(trace(vm\c));
    end
end
end