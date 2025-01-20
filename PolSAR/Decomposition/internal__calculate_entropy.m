function h = internal__calculate_entropy(x)

h = 0;
len = numel(x);
for i = 1:len
    if x ~= 0
        h = h - x(i) * log(x(i)) / log(len);
    end
end

end