function [H,alpha,A] = CloudePottier(T3)

arguments
    T3 (:,:,3,3)
end

[height,width,~,~] = size(T3);
T3 = parallel.pool.Constant(shiftdim(T3, 2));

H = zeros(height, width);
alpha = zeros(height, width);
A = zeros(height, width);
parfor j = 1:width
    for i = 1:height
        t = squeeze(T3.Value(:,:,i,j));
        [v,d] = eig(t);
        d = diag(abs(d))'
        p = d / sum(d);
        H(i,j) = -sum(p .* log(p) / log(3));
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        l3 = min(d);
        l2 = sum(d) - max(d) - l3;
        A(i,j) = (l2 - l3) / (l2 + l3);
    end
end

end