function [H,alpha,A] = internal__CloudePottier_matlab(M3)

height = M3.Height;
width = M3.Width;

H = zeros(height, width, M3.Dtype);
alpha = zeros(height, width, M3.Dtype);
A = zeros(height, width, M3.Dtype);

M3 = parallel.pool.Constant(M3);
parfor j = 1:width
    for i = 1:height
        t = M3.Value.MatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = internal__calculate_entropy(p);
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        l3 = min(d);
        l2 = sum(d) - max(d) - l3;
        A(i,j) = (l2 - l3) / (l2 + l3);
    end
end

end