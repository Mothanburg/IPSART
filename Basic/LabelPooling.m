% LabelPooling - 根据区域索引标签的池化
% 图像的区域索引标签必须是连续的 1..N
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function result = LabelPooling(image, label, method)

arguments
    image (:,:)
    label (:,:)
    method = @mean % 默认方法为均值池化（Average Pooling）
end

vals = accumarray(label(:), image(:), [max(label(:)) 1], method);
result = vals(label);

end