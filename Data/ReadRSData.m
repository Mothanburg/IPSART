% ReadRSData - 兼容ENVI格式的数据读取器
% fileExtensionPolicy参数用于指定数据头文件的命名方式
%   "SUB" 表示数据文件"a.x"对应的头文件名称为"a.hdr"
%   "APPEND" 表示数据文件"b.y"对应的头文件名称为"b.y.hdr"
% 如果数据文件没有后缀名，请使用"APPEND"
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function data = ReadRSData(filePath, options)

arguments
    filePath string
    options.fileExtensionPolicy string = "Sub"
end

[prefix,file_name,file_ext] = fileparts(filePath);

switch upper(options.fileExtensionPolicy)
    case "SUB"
        header_file = fullfile(prefix, file_name + ".hdr");
    case "APPEND"
        header_file = fullfile(prefix, file_name + file_ext + ".hdr");
    otherwise
        error("Unknown option of fileExtensionPolicy: '%s'", ...
            options.fileExtensionPolicy);
end

if ~exist(filePath, "file")
    error("Cannot find the data file '%s'", file_name + file_ext);
elseif ~exist(header_file, "file")
    error("Cannot find the header file of '%s'", file_name + file_ext);
end

header = fn_parse_envi_header(header_file);

if startsWith(header.data_type, "complex")

    inter_datatype = replace(header.data_type, "complex", "float");
    raw = multibandread(...
        filePath, ...
        [header.lines header.samples * 2 header.bands], ...
        "*" + inter_datatype, ...
        header.header_offset, ...
        header.interleave, ...
        header.byte_order ...
        );
    bands = cell(1, header.bands);
    for i = 1:header.bands
        band_real = raw(:,1:2:end,i);
        band_imag = raw(:,2:2:end,i);
        bands{i} = squeeze(band_real + 1i * band_imag);
    end
    data = cat(3, bands{:});

else

    data = multibandread(...
        filePath, ...
        [header.lines header.samples header.bands], ...
        "*" + header.data_type, ...
        header.header_offset, ...
        header.interleave, ...
        header.byte_order ...
        );

end

end


% Parsing ENVI header file
% TODO: support more properties in ENVI header file
% NOTE: This function may cause error because of nonstandard syntax in the file
function headerInfo = fn_parse_envi_header(headerfile)

lines = readlines(headerfile, "WhitespaceRule", "trim", "EmptyLineRule", "skip");
assert(strcmp(lines(1), "ENVI"), "The first line of .hdr file must be 'ENVI'");

headerInfo = struct();

n_line = 2;
while true
    line = char(lines(n_line));

    % Parse property name
    [token,n_char] = get_next_token(line, 1);
    if token.type ~= "WORD"
        error("The property name must start with a letter.");
    end

    property_name = token.value;
    % ENVI headers may contain property names with multiple words (e.g., "map info").
    % These are represented as separate WORD tokens and joined with underscores.
    while true
        [token,n_char] = get_next_token(line, n_char);
        if token.type ~= "WORD"
            break
        end
        property_name = property_name + "_" + token.value;
    end

    % Next token must be 'ASSIGN'
    if token.type ~= "ASSIGN"
        error("Missing '=' in line %d.", n_line);
    end

    % Parse property value
    [token,n_char] = get_next_token(line, n_char);
    if token.type == "STRING" || token.type == "WORD"
        property_value = token.value;
        while true
            [token,n_char] = get_next_token(line, n_char);
            if token.type == "NULL"
                break
            end
            property_value = property_value + " " + token.value;
        end
    elseif token.type == "L_BRACE"
        brace_pos = strfind(line, '}');
        if ~isempty(brace_pos) && brace_pos(1) >= n_char
            content = line(n_char:brace_pos(1)-1);
        else
            content = line(n_char:end);
            n_line = n_line + 1;
            while n_line <= length(lines)
                next_line = char(lines(n_line));
                brace_pos = strfind(next_line, '}');
                if ~isempty(brace_pos)
                    content = [content '\n' next_line(1:brace_pos(1)-1)];
                    break
                else
                    content = [content '\n' next_line];
                    n_line = n_line + 1;
                end
            end
            if n_line > length(lines)
                error("Unclosed '}' in .hdr file");
            end
        end
        property_value = string(content);
    else
        error("Invalid property value in line %d.", n_line);
    end

    headerInfo.(property_name) = fn_lookup_valuetype(property_name, property_value);
    n_line = n_line + 1;
    if n_line > length(lines)
        break
    end
end


end

function [token,ptr] = get_next_token(text, start_ptr)

idx = start_ptr;
while idx <= length(text) && isspace(text(idx))
    idx = idx + 1;
end

% valid type: NULL (no more token) WORD (starts with letter), ASSIGN ('='), L_BRACE ('{'),
% STRING (any non-space ascii character that is not a WORD)

if idx > length(text)
    token.type = "NULL";
    ptr = idx;
elseif text(idx) == '='
    token.type = "ASSIGN";
    ptr = idx + 1;
elseif text(idx) == '{'
    token.type = "L_BRACE";
    ptr = idx + 1;
elseif isletter(text(idx))
    content = '';
    while idx <= length(text) && ~isspace(text(idx))
        content = [content text(idx)];
        idx = idx + 1;
    end
    token.type = "WORD";
    token.value = string(content);
    ptr = idx;
else
    content = '';
    while idx <= length(text) && ~isspace(text(idx))
        content = [content text(idx)];
        idx = idx + 1;
    end
    token.type = "STRING";
    token.value = string(content);
    ptr = idx;
end


end

% Convert the property name to property value
function value = fn_lookup_valuetype(property_name, value_string)

switch property_name
    case "bands"  % band num
        value = str2double(value_string);
    case "byte_order"  % byte order, little edian: 0, big edain: 1
        if ismember(value_string, ["ieee-le", "ieee-be"])
            value = value_string;
        else
            edians = ["ieee-le" "ieee-be"];
            value = edians(str2double(value_string) + 1);
        end
    case "data_type"
        type = str2double(value_string);
        switch type
            case 1
                value = "uint8";
            case 2
                value = "int16";
            case 3
                value = "int32";
            case 4
                value = "float32";
            case 5
                value = "float64";
            case 6
                value = "complex32";
            case 9
                value = "complex64";
            case 12
                value = "uint16";
            case 13
                value = "uint32";
            case 14
                value = "int64";
            case 15
                value = "uint64";
            otherwise
                error("Unknown data type");
        end
    case "header_offset"
        value = str2double(value_string);
    case "interleave"
        validatestring(value_string, ["bsq", "bil", "bip"], ...
            'fn_lookup_valuetype', 'interleave');
        value = value_string;
    case "samples"
        value = str2double(value_string);
    case "lines"
        value = str2double(value_string);
    otherwise
        value = value_string;
end

end

