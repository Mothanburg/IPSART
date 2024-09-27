function [init_classes, num_classes] = gen_init_cls_Freeman(Ps, Pd, Pv, C, numClasses)

[height,width] = size(Ps);
% 初始化30个类别
[~,scatter_classes] = max(abs(cat(3, Ps, Pd, Pv)), [], 3);
surface = scatter_classes==1;
double = scatter_classes==2;
volume = scatter_classes==3;
s_borders = divide_vector(sort(Ps(surface)), 10);
d_borders = divide_vector(sort(Pd(double)), 10);
v_borders = divide_vector(sort(Pv(volume)), 10);
init_classes = zeros(height, width);
for m=1:height
    for n=1:width
        switch scatter_classes(m,n)
            case 1
                init_classes(m,n) = 10;
                for idx=1:9
                    if Ps(m,n)<=s_borders(idx)
                        init_classes(m,n) = idx;
                        break;
                    end
                end
            case 2
                init_classes(m,n) = 20;
                for idx=1:9
                    if Pd(m,n)<=d_borders(idx)
                        init_classes(m,n) = 10+idx;
                        break;
                    end
                end
            case 3
                init_classes(m,n) = 30;
                for idx=1:9
                    if Pv(m,n)<=v_borders(idx)
                        init_classes(m,n) = 20+idx;
                        break;
                    end
                end
        end
    end
end
% 合并类别
num_classes = numClasses;
merged_classes = init_classes;
s_classes = 1:10;
d_classes = 11:20;
v_classes = 21:30;
for cnt=30:-1:(num_classes+1)
    c = calc_centers(cnt, C, merged_classes);
    distances = Inf(cnt);
    for i=1:(cnt-1)
        for j=(i+1):cnt
            d = (log(real(det(c{i})))+log(real(det(c{j})))+...
                trace(c{i}\c{j}+c{j}\c{i}))/2;
            distances(i,j) = d;
        end
    end
    % 确保合并的类别属于一个散射体制，同时总数不会过多
    while 1
        [~,I] = min(distances, [], "all");
        [p,q] = ind2sub(size(distances), I); % 计算出类别p和类别q
        % p和q的散射机制相同
        is_surface = ismember(p, s_classes) && ismember(q, s_classes);
        is_double = ismember(p, d_classes) && ismember(q, d_classes);
        is_volume = ismember(p, v_classes) && ismember(q, v_classes);
        if is_surface || is_double || is_volume
            % p和q合并后数量不过多
            if sum(merged_classes==p, "all")+...
                    sum(merged_classes==q, "all")<=2*height*width/(cnt-1)
                break;
            end
        end
        distances(p,q) = inf;
    end
    % 合并类别
    for m=1:height
        for n=1:width
            if merged_classes(m,n)==q
                merged_classes(m,n) = p;
            elseif merged_classes(m,n)>q
                merged_classes(m,n) = merged_classes(m,n)-1;
            end
        end
    end
    % 更新三种散射机制的类别表
    if is_surface
        s_classes(s_classes>=q) = s_classes(s_classes>=q)-1;
        s_classes = sort(unique(s_classes));
        d_classes = d_classes-1;
        v_classes = v_classes-1;
    elseif is_double
        d_classes(d_classes>=q) = d_classes(d_classes>=q)-1;
        d_classes = sort(unique(d_classes));
        v_classes = v_classes-1;
    elseif is_volume
        v_classes(v_classes>=q) = v_classes(v_classes>=q)-1;
        v_classes = sort(unique(v_classes));
    end
end
init_classes = merged_classes;

end


% 在排好序的向量中等间隔地取N个值
function borders = divide_vector(V, N)

borders = zeros(1,N);
len = length(V);
for i=1:N
    borders(i) = V(round(len*i/N));
end

end
