function [H,alpha,A] = internal__DualPolCloudePottier_matlab(M2)

height = M2.Height;
width = M2.Width;

H = zeros(height, width, M2.Dtype);
alpha = zeros(height, width, M2.Dtype);
A = zeros(height, width, M2.Dtype);

M2 = parallel.pool.Constant(M2);
parfor j = 1:width
    for i = 1:height
        t = M2.Value.MatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = internal__calculate_entropy(p);
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        A(i,j) = abs(d(1) - d(2)) / sum(d);
    end
end

end
