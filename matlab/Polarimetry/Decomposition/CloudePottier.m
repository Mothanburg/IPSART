function [H,alpha,A] = CloudePottier(T3)

arguments
    T3 PolT3
end

% global GARS_CONFIG

[H,alpha,A] = CloudePottierL0(T3);

end


function [H,alpha,A] = CloudePottierL0(T3)

height = T3.Height;
width = T3.Width;

H = zeros(height, width, T3.Dtype);
alpha = zeros(height, width, T3.Dtype);
A = zeros(height, width, T3.Dtype);

T3 = parallel.pool.Constant(T3);
parfor j = 1:width
    for i = 1:height
        t = T3.Value.getMatAt(i, j);
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