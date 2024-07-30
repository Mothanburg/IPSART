function composed = ComposeOp(funcs)

arguments (Repeating)
    funcs function_handle
end

composed = @(x) seq_call(x, funcs);

end

function result = seq_call(x, funcs)

result = x;
for fn = funcs
    result = feval(fn{:}, result);
end

end