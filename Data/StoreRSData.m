% Store data as ENVI data format
% File name of stored data and its header file:
%    "prefix/dataName" and "prefix/dataName.hdr"
function data = StoreRSData(data, prefix, dataName, options)

arguments
    data (:,:,:)
    prefix string
    dataName string = ""
    options.fileExt string = ""
    options.useLowPrecision logical = false
end

if ~exist(prefix, "dir")
    mkdir(prefix);
end

if ~endsWith(prefix, [filesep, "/"])
    prefix = strcat(prefix, filesep);
end

if dataName == ""
    dataName = inputname(1);
end

filename = strcat(prefix, dataName, options.fileExt);
if exist(filename, "file")
    warning("The file %s alread exists, overwriting it.", filename);
    delete(filename);
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
        filename, ...
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
        filename, ...
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
fid = fopen(strcat(prefix, dataName, ".hdr"), "w");
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

