% ComposeOp - 将多个单参数函数（一元谓词）按顺序组合成一个
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function composed = ComposeOp(unaries)

arguments (Repeating)
    unaries function_handle
end

composed = @(x) seq_call(x, unaries);

end

function result = seq_call(x, funcs)

result = x;
for fn = funcs
    result = feval(fn{:}, result);
end

end