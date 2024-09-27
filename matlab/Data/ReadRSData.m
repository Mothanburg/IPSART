% Read Data with ENVI header file
% Corresponding file names of data and header file: 
%     "prefix/dataName.fileExt" and "prefix/dataName.hdr"
function data = ReadRSData(prefix, dataName, fileExt)

arguments
    prefix string
    dataName string
    fileExt string = ""
end

if ~endsWith(prefix, [filesep, "/"])
    prefix = strcat(prefix, filesep);
end
file_path = strcat(prefix, dataName, fileExt);
header_path = strcat(prefix, dataName, ".hdr");
if ~exist(header_path, "file") || ~exist(file_path, "file")
    error("Cannot find data or its .hdr file");
end

header = fn_parse_envi_header(header_path);

if startsWith(header.data_type, "complex")

    inter_datatype = replace(header.data_type, "complex", "float");
    raw = multibandread(...
        file_path, ...
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
        file_path, ...
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
    if token.type == "VALUE" || token.type == "WORD"
        property_value = token.value;
        while true
            [token,n_char] = get_next_token(line, n_char);
            if token.type == "NULL"
                break
            end
            property_value = property_value + " " + token.value;
        end
    elseif token.type == "L_BRACE"
        if line(end) == '}'
            content = line(n_char:end-1);
        else
            content = line(n_char:end);
            n_line = n_line + 1;

            while n_line <= length(lines) && ~endsWith(lines(n_line), "}")
                content = ['\n' content char(lines(n_line))];
                n_line = n_line + 1;
            end

            if n_line > length(lines)
                error("Unclosed '}' in .hdr file");
            end

            content = [content char(lines(n_line))];
        end
        property_value = string(content(1:end-1));
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
% VALUE (any non-space ascii character)

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
    token.type = "VALUE";
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
        edians = ["ieee-le" "ieee-be"];
        value = edians(str2double(value_string) + 1);
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
        assert("bsq" == value_string || "bil" == value_string || "bip" == value_string)
        value = value_string;
    case "samples"
        value = str2double(value_string);
    case "lines"
        value = str2double(value_string);
    otherwise
        value = value_string;
end

end

