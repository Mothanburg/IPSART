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