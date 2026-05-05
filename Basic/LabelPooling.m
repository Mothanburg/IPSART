% 根据索引标签划分的图像区域进行池化操作
% label必须是连续的 1..N
% 默认为平均值池化
function result = LabelPooling(image, label, method)

arguments
    image (:,:)
    label (:,:)
    method = @mean 
end

vals = accumarray(label(:), image(:), [max(label(:)) 1], method);
result = vals(label);

end