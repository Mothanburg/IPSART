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
        t = T3.Value(:,:,i,j);
        [u,l] = eig(t);
        l = sort(abs(diag(l)'), 'descend');
        p = l / trace(t);
        H(i,j) = -sum(p .* log(p)) / log(3);
        alpha(i,j) = 180 * sum(p .* acos(abs(u(1,:)))) / pi;
        A(i,j) = (l(2) - l(3)) / (l(2) + l(3));
    end
end

end