% 全极化H/a/A的初始类别
function [init_classes, num_classes] = gen_init_cls_H_a_A(H, a, A)

[init_classes,~] = gen_init_cls_H_a(H, a);
num_classes = 16;
init_classes = init_classes+(8*(A>0.5));

end