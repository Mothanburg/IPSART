function result = PolTransform(from, to, varargin)

arguments
    from string
    to string
end

arguments (Repeating)
    varargin
end

switch from
    case {"", "Gen"}
        switch to
            case "S2"
                transformer = @fn_gen_S2;
            case "C3"
                transformer = @fn_gen_C3;
            case "T3"
                transformer = @fn_gen_T3;
            otherwise
                error("PolTransform:invalidOutputFormat", "错误的输出格式");
        end
    case "S2"
        switch to
            case "C3"
                transformer = @fn_S2_to_C3;
            case "T3"
                transformer = @fn_S2_to_T3;
            otherwise
                error("PolTransform:invalidOutputFormat", "错误的输出格式");
        end
    case "C3"
        switch to
            case "T3"
                transformer = @fn_C3_to_T3;
            otherwise
                error("PolTransform:invalidOutputFormat", "错误的输出格式");
        end
    case "T3"
        switch to
            case "C3"
                transformer = @fn_T3_to_C3;
            otherwise
                error("PolTransform:invalidOutputFormat", "错误的输出格式");
        end
        % case "C2.HHVV"
        %     transformer = ["T2.HHVV"];
        % case "C2.HHHV"
        %     transformer = ["T2.HHHV"];
        % case "C2.VVVH"
        %     transformer = ["T2.VVVH"];
        % case "CP.PI4"
        %     transformer = [];
        % case "CP.CTLR"
        %     transformer = [];
        % case "CP.DCP"
        %     transformer = [];
    otherwise
        error("PolTransform:invalidInputFormat", "错误的输入格式");
end

try
    [path,~,~] = fileparts(mfilename("fullpath"));
    addpath(strcat(path, filesep, "impl"));
    result = transformer(varargin{:});
catch cause
    rmpath(strcat(path, filesep, "impl"));
    err = MException("PolTransform:cannotTransform", "无法完成极化变换");
    err.addCause(cause);
    throw(err);
end
rmpath(strcat(path, filesep, "impl"));

end

