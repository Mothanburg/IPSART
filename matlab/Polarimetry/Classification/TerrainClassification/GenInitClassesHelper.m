% 生成初始类别
function [initClasses, numClasses] = GenInitClassesHelper(method, input)

arguments
    method string
    input.H (:,:) double = []
    input.a (:,:) double = []
    input.A (:,:) double = []
    input.PD (:,:) double = []
    input.C (:,:) double = []
    input.Ps (:,:) double = []
    input.Pd (:,:) double = []
    input.Pv (:,:) double = []
    input.numClasses = 8
end

switch method
    case 'H/a'
        [initClasses, numClasses] = gen_init_cls_H_a(input.H, input.a);
    case 'H/a/A'
        [initClasses, numClasses] = gen_init_cls_H_a_A(input.H, input.a, input.A);
    case 'Dual H/a'
        [initClasses, numClasses] = gen_init_cls_Dual_H_a(input.H, input.a);
    case 'Freeman'
        [initClasses, numClasses] = ...
            gen_init_cls_Freeman(input.Ps, input.Pd, input.Pv, input.C, input.numClasses);
    otherwise
        error('不支持的操作');
end

end