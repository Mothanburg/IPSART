% internal__UnsupervisedG0Clustering_matlab - Internal MATLAB implementation of G0 clustering
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function clustered = internal__UnsupervisedG0Clustering_matlab( ...
    C, lookNum, initClasses, numClasses, threshold, maxIter)

height = C.Height;
width = C.Width;

transfer_rate = 1;
iter_input = initClasses;

iters = 1;
len = height * width;
while transfer_rate > threshold && iters <= maxIter
    % Calculate center cov matrix and esimate alpha parameter
    vms = zeros(C.Dim, C.Dim, numClasses);
    alphas = zeros(C.Dim, numClasses);
    for i = 1:C.Dim
        for j = i:C.Dim
            page = C.PageAt(i,j);
            if i == j
                for cls = 1:numClasses
                    m1 = mean(page(iter_input == cls));
                    m2 = mean(page(iter_input == cls).^2);
                    alphas(i,cls) = m1^2 * (lookNum + 1) / (m2 * lookNum - m1^2 * (lookNum + 1));
                    vms(i,j,cls) = mean(page(iter_input == cls));
                end
            else
                for cls = 1:numClasses
                    vms(i,j,cls) = mean(page(iter_input == cls));
                    vms(j,i,cls) = conj(vms(i,j,cls));
                end
            end
        end
    end
    alphas = mean(alphas, 1);

    % update label
    iter_output = zeros(height, width);
    parfor idx = 1:len
        [i,j] = ind2sub([height, width], idx);
        c = C.MatAt(i, j);
        distances = calc_distances(c, vms, alphas, lookNum);
        [~,new_class] = min(distances);
        iter_output(idx) = new_class;
    end
    transfer_rate = sum(iter_output~=iter_input, 'all') / len;

    iter_input = iter_output;
    iters = iters + 1;
end
clustered = iter_output;

end


function ds = calc_distances(c, vms, alphas, n)
num_classes = size(vms, 3);
ds = zeros(1, num_classes);
for m = 1:num_classes
    vm = squeeze(vms(:,:,m));
    if any(isnan(vm(:)))
        ds(m) = inf;
    else
        q = size(c, 1);
        a = alphas(m);
        trace_term = real(trace(vm \ c));
        ds(m) = (n + 1) * log(real(det(vm)))...
            + a * log(-a - 1) ...
            - log(gamma(q * n - a)) ...
            + log(gamma(-a)) ...
            + n * trace_term ...
            - (a - q * n) * log(n * trace_term - a - 1);
    end
end
end