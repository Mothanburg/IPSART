function [H,a,A] = CloudePottier4DP(T2)

arguments
    T2 (:,:,2,2)
end

[height,width,~,~] = size(T2);
H = zeros(height, width);
a = zeros(height, width);
A = zeros(height, width);

T2 = parallel.pool.Constant(shiftdim(T2, 2));
parfor j = 1:width
    for i = 1:height
        t = T2.Value(:,:,i,j);
        [u,l] = eig(t);
        l = abs(diag(l));
        p = l / trace(t);
        H(i,j) = -sum(p .* log2(p));
        a(i,j) = 180 * sum(p .* acos(abs(u(1,:)))) / pi;
        A(i,j) = abs(l(1) - l(2)) / trace(t);
    end
end

end
