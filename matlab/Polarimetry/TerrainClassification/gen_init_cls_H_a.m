% 全极化H/a平面的初始类别
function [init_classes, num_classes] = gen_init_cls_H_a(H, a)

[height,width] = size(H);
num_classes = 8;
init_classes = zeros(height, width);
for i = 1:height
    for j = 1:width
        if H(i,j) <= 0.5
            if a(i,j) >= 47.5
                init_classes(i,j) = 1;
            elseif a(i,j) >= 42.5
                init_classes(i,j) = 2;
            else
                init_classes(i,j) = 3;
            end
        elseif H(i,j) <= 0.9
            if a(i,j) >= 50
                init_classes(i,j) = 4;
            elseif a(i,j) >= 40
                init_classes(i,j) = 5;
            else
                init_classes(i,j) = 6;
            end
        else
            if a(i,j) >= 55
                init_classes(i,j) = 7;
            elseif a(i,j) >= 40
                init_classes(i,j) = 8;
            else
                init_classes(i,j) = 9; % 不存在
            end
        end
    end
end

end