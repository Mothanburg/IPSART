function InitGaRSPackages(packageList)

arguments
    packageList (1,:) string = ["基础工具" "极化SAR" "图像处理"]
end

[path,~,~] = fileparts(mfilename("fullpath"));
for name = packageList
    fn_add_package(path, name);
end

assignin("base", "L1PSH", @(x) HistStretch(x, "Linear Percent", 0, 99));
assignin("base", "L2PSH", @(x) HistStretch(x, "Linear Percent", 0, 98));
assignin("base", "L5PSH", @(x) HistStretch(x, "Linear Percent", 0, 95));

end


function fn_add_package(path, name)

try
    package_path = strcat(path, filesep, name);
    fn_add_functions(package_path);
catch cause
    rehash("path");
    err = MException("GaRS:ImportPackages:cannotImportPackage", "未能导入函数包：%s", name);
    err.addCause(cause);
    throw(err);
end

end


function fn_add_functions(package_path)

if ~exist(package_path, "dir")
    error("文件夹""%s""不存在");
end

file_list = dir(package_path);

for i = 1:length(file_list)
    if file_list(i).isdir && ~strcmp(file_list(i).name, ".") && ~strcmp(file_list(i).name, "..")
        addpath(strcat(file_list(i).folder, filesep, file_list(i).name));
    end
end

end