function composed = ComposeOp(funcs)

arguments (Repeating)
    funcs function_handle
end

composed = @(x) seq_call(x, funcs);

end

function result = seq_call(x, funcs)

if isempty(funcs)
    result = x;
else
    result = seq_call(feval(funcs{1}, x), funcs(2:end));
end

end