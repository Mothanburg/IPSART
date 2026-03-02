function [H,alpha,A] = internal__CloudePottier_matlab(M)

height = M.Height;
width = M.Width;

H = zeros(height, width, M.Dtype);
alpha = zeros(height, width, M.Dtype);
A = zeros(height, width, M.Dtype);

parfor j = 1:width
    for i = 1:height
        t = M.MatAt(i, j);
        [v,d] = eig(t);
        d = diag(abs(d))';
        p = d / sum(d);
        H(i,j) = calculate_entropy(p);
        alpha(i,j) = sum(p .* acosd(abs(v(1,:))));
        if M.Dim == 3
            l3 = min(d);
            l2 = sum(d) - max(d) - l3;
            A(i,j) = (l2 - l3) / (l2 + l3);
        else
            A(i,j) = p(1) - p(2);
        end
    end
end

end


function h = calculate_entropy(x)
h = 0;
len = numel(x);
for i = 1:len
    if x ~= 0
        h = h - x(i) * log(x(i)) / log(len);
    end
end
end