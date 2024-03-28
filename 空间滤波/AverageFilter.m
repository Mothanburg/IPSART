% 均值滤波，不改变图像尺寸的多视处理
function result = AverageFilter(image, rowLook, colLook)

arguments
    image double
    rowLook double
    colLook double
end

sz = size(image);
if length(sz) < 2
    error("输入值至少必须是二维图像");
end

h = fspecial("average", [rowLook colLook]);

if length(sz) > 2
    page_cnt = prod(sz(3:end));
    result = zeros(sz);
else
    page_cnt = 1;
    image = reshape(image, [sz 1]);
    result = zeros([sz 1]);
end

for idx = 1:page_cnt
    result(:,:,idx) = imfilter(squeeze(image(:,:,idx)), h);
end

result = squeeze(result);

end