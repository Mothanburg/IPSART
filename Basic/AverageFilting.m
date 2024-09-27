% 均值滤波，不改变图像尺寸的多视处理
function result = AverageFilting(image, rowLook, colLook)

arguments
    image (:,:) double
    rowLook double
    colLook double
end

h = fspecial("average", [rowLook colLook]);
result = imfilter(image, h);

end