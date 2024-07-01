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
    error("ReadRSData:cannotFindFile", "找不到数据文件或其头文件");
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


% 解析ENVI头文件，允许加入自定义属性
% 语法要求：属性名 + <任意数量空格/制表符> + "="符号 + <任意数量空格/制表符> + 属性值 + 换行符
% 属性名要求：全字母，单词间仅允许单个空格，解析后空格会被替换为"_"
% 属性值要求：一般属性不允许包含"="、"；"，且必须在单行内写完；如要包含特殊字符或换行，请用"{}"将属性值括起来
% TODO：完整支持标准ENVI头文件
function headerInfo = fn_parse_envi_header(headerfile)

lines = readlines(headerfile, "WhitespaceRule", "trim", "EmptyLineRule", "skip");
assert(strcmp(lines(1), "ENVI"), "ReadRSData:invalidHeaderFormat", "头文件第一行必须是'ENVI'");

raw_text = char(join(lines(2:end), ';') + ';');

headerInfo = struct();

% 状态机
% 0：初始状态
% 1：解析属性名状态
% 2：解析属性值状态
state = 0;
i = 1;
while true

    switch state
        case 0
            assert(isletter(raw_text(i)), "ReadRSData:invalidPropertyName", "属性名开头必须是字母");
            state = 1;

        case 1  % 解析属性名
            property_name = '';
            while true
                if isletter(raw_text(i))
                    property_name = [property_name raw_text(i)];
                    i = i + 1;
                elseif raw_text(i) == ' ' && isletter(raw_text(i + 1))
                    property_name = [property_name '_' raw_text(i + 1)];
                    i = i + 2;
                else
                    break;
                end
            end

            % 跳过属性名后的空白
            while isspace(raw_text(i))
                i = i + 1;
            end

            assert(raw_text(i) == '=', "ReadRSData:headerSyntaxError", "属性名后同一行内的下一个token必须是'='符号")
            i = i + 1;

            % 跳过Value前的空白
            while isspace(raw_text(i))
                i = i + 1;
            end

            state = 2;

        case 2  % 解析属性值
            value_string = '';
            % 暂时不解析用"{ }"包含的内容，直接将这些内容转换为字符串
            if raw_text(i) == '{'
                i = i + 1;
                while raw_text(i) ~= '}'
                    if raw_text(i) == ';'
                        value_string = [value_string '\n'];
                    else
                        value_string = [value_string raw_text(i)];
                    end
                    i = i + 1;
                end
                i = i + 1;
                headerInfo.(property_name) = value_string;
                % 解析普通的属性值
            else
                while raw_text(i) ~= ';'
                    value_string = [value_string raw_text(i)];
                    i = i + 1;
                end

                value = fn_lookup_valuetype(property_name, value_string);
                headerInfo.(property_name) = value;
            end

            if i == length(raw_text)
                break;
            else
                i = i + 1;
                state = 1;
            end
    end

end

end

% 根据ENVI头文件属性名，将属性值转换为真实类型
% NOTE：现只支持必要的属性，其它属性将保留成字符串
% TODO：支持全部属性
function value = fn_lookup_valuetype(property_name, value_string)

switch property_name
    case "bands"  % 波段数
        value = str2double(value_string);
    case "byte_order"  % 字节序，小端：0，大端：1
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
                error("ReadRSData:invalidDataType", "错误的数据类型");
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

