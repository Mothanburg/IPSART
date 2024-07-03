function wrapped = AsUnary(func, position, vargin)

arguments
    func function_handle
    position {mustBeInteger, mustBeGreaterThanOrEqual(position, 1)}
end

arguments (Repeating)
    vargin
end

wrapped = @(x) func(vargin{1:(position - 1)}, x, vargin{position:end});

end