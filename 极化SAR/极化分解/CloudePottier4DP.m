function [H,alpha,A] = CloudePottier4DP(T2)

arguments
    T2 (:,:,2,2)
end

[height,width,~,~] = size(T2);
H = zeros(height, width);
alpha = zeros(height, width);
A = zeros(height, width);

T2 = parallel.pool.Constant(shiftdim(T2, 2));
parfor j = 1:width
    for i = 1:height
        t = squeeze(T2.Value(:,:,i,j));
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = -sum(p .* log2(p));
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        A(i,j) = abs(d(1) - d(2)) / sum(d);
    end
end

end
