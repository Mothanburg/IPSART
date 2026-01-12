% Store data as ENVI data format
% Parameters
%    data: array with up to 3 dimensions
%    filePath: full path of the data
%    fileExtensionPolicy: Determine the file name of the generated header file
%       "Sub" [default]: substitution the data file's extension into .hdr
%       "Append": append .hdr to the data file's name 
%    useLowPrecision: whether cast double/int64 to single/int32 when storing 
function StoreRSData(data, filePath, options)

arguments
    data (:,:,:)
    filePath string
    options.fileExtensionPolicy string = "Sub"
    options.useLowPrecision logical = false
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

if exist(filePath, "file")
    warning("The file %s exists, overwriting it.", filePath);
    delete(filePath);
elseif ~exist(prefix, "dir") && prefix ~= ""
    mkdir(prefix);
end

[lines,samples,bands] = size(data);
if ~isreal(data)
    if class(data) == "double" && ~options.useLowPrecision
        data_type = "complex64";
    else
        % ENVI doesn't support non-float complex number
        data_type = "complex32";
    end
    
    raw_bands = zeros(lines, samples * 2, bands);
    raw_bands(:,1:2:end,:) = real(data);
    raw_bands(:,2:2:end,:) = imag(data);
    
    multibandwrite(...
        cast(raw_bands, class(data)), ...
        filePath, ...
        "bsq", ...
        "machfmt", "ieee-le" ...
        );
else
    data_type = class(data);
    
    if options.useLowPrecision
        switch data_type
            case {"uint32", "uint64"}
                data_type = "uint16";
                data = cast(data, data_type);
            case {"int32", "int64"}
                data_type = "int16";
                data = cast(data, data_type);
            case "double"
                data_type = "single";
                data = cast(data, data_type);
        end
    end
    
    multibandwrite(...
        data, ...
        filePath, ...
        "bsq", ...
        "machfmt", "ieee-le" ...
        );
end

switch data_type
    case "uint8"
        data_type_v = 1;
    case "int16"
        data_type_v = 2;
    case "int32"
        data_type_v = 3;
    case "single"
        data_type_v = 4;
    case "double"
        data_type_v = 5;
    case "complex32"
        data_type_v = 6;
    case "complex64"
        data_type_v = 9;
    case "uint16"
        data_type_v = 12;
    case "uint32"
        data_type_v = 13;
    case "int64"
        data_type_v = 14;
    case "uint64"
        data_type_v = 15;
    otherwise
        error('Unknown error')
end

% Write header file
fid = fopen(header_file, "w");
fprintf(fid, "ENVI\n");
fprintf(fid, "bands = %d\n", bands);
fprintf(fid, "samples = %d\n", samples);
fprintf(fid, "lines = %d\n", lines);
fprintf(fid, "header offset = 0\n");
fprintf(fid, "interleave = bsq\n");
fprintf(fid, "byte order = 0\n");
fprintf(fid, "file type = ENVI Standard\n");
fprintf(fid, "data type = %d\n", data_type_v);
fclose(fid);

end

