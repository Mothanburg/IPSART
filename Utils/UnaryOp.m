% UnaryOp - 将函数封装为单参数函数（一元谓词）
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function unary = UnaryOp(func, position, vargin)

arguments
    func function_handle
    position {mustBeInteger, mustBeGreaterThanOrEqual(position, 1)}
end

arguments (Repeating)
    vargin
end

unary = @(x) func(vargin{1:(position - 1)}, x, vargin{position:end});

end