% 全极化 SAR 图像的二分分区树（Binary Partition Tree，BPT）构建
classdef PolBPT

    properties
        Us
        Vs
        Ps
        Costs
        Phis
        Root
        ImageSize
    end

    properties (Dependent)
        NumLeaves
    end

    methods
        % 构造函数，使用极化协方差矩阵构建 BPT 树
        function obj = PolBPT(M)
            arguments
                M PolMat
            end

            obj.ImageSize = [M.Height, M.Width];
            n = M.Height * M.Width;

            % 返回 mergers 矩阵: [u, v, p, cost, phi]
            % mergers从第一列到最后一列，记录了每一次合并的过程
            global IPSARTMexHost;
            if isempty(IPSARTMexHost)
                IPSARTMexHost = mexhost();
            end
            if M.Dim == 2
                [obj.Us, obj.Vs, obj.Ps, obj.Costs, obj.Phis] = ...
                    IPSARTMexHost.feval("internal__mex_bridge", "PolBPT::Build", ...
                    M.m11, M.m22, M.m12_r, M.m12_i);
            elseif M.Dim == 3
                [obj.Us, obj.Vs, obj.Ps, obj.Costs, obj.Phis] = ...
                    IPSARTMexHost.feval("internal__mex_bridge", "PolBPT::Build", ...
                    M.m11, M.m22, M.m33, M.m12_r, M.m13_r, M.m23_r, M.m12_i, ...
                    M.m13_i, M.m23_i);
            elseif M.Dim == 4
                [obj.Us, obj.Vs, obj.Ps, obj.Costs, obj.Phis] = ...
                    IPSARTMexHost.feval("internal__mex_bridge", "PolBPT::Build", ...
                    M.m11, M.m22, M.m33, M.m44, M.m12_r, M.m13_r, M.m14_r, ...
                    M.m23_r, M.m24_r, M.m34_r, M.m12_i, M.m13_i, M.m14_i, ...
                    M.m23_i, M.m24_i, M.m34_i);
            else
                error("Unknown error: undefined dim of PolMat object.")
            end

            % 根节点是最后一次合并产生的新节点
            if ~isempty(obj.Ps)
                obj.Root = obj.Ps(end);
            else
                obj.Root = n;
            end
        end

        function value = get.NumLeaves(obj)
            value = prod(obj.ImageSize);
        end

        % 根据想要的区域数量对 BPT 进行剪枝，并输出图像标签
        function label = PruneByNumRegions(obj, target_num_regions)

            % 在 BPT 的构建过程中，我们进行了 num_leaves - 1 次合并
            % 如果想要保留 K 个区域，我们需要保留前 (num_leaves - K) 次合并
            % 也就是撤销最后的 (K - 1) 次合并
            num_steps_to_keep = obj.NumLeaves - target_num_regions;

            if num_steps_to_keep < 0
                warning('目标区域数大于初始超像素数，返回原始标签');
                label = reshape(1:obj.NumLeaves, obj.ImageSize);
                return;
            end

            % 剪枝后，BPT 树变成了一个森林，需要自底向上合并到每棵树的树根
            % 对此为了性能，我们使用路径压缩算法找到每个像素（叶子）的最终树根
            map_idx = 1:obj.Root;
            us = obj.Us(1:num_steps_to_keep);
            vs = obj.Vs(1:num_steps_to_keep);
            ps = obj.Ps(1:num_steps_to_keep);
            map_idx(us) = ps;
            map_idx(vs) = ps;

            global IPSARTMexHost;
            if isempty(IPSARTMexHost)
                IPSARTMexHost = mexhost();
            end
            final_idx = IPSARTMexHost.feval("internal__mex_bridge", ...
                "PolBPT::FindRoot", cast(map_idx, 'int32'));

            % 前 1..num_leaves 个结点是每个像素的标签
            pixel_idx = final_idx(1:obj.NumLeaves);

            [~, ~, label] = unique(pixel_idx);
            label = reshape(label, obj.ImageSize);
        end

        % 根据图像同质度阈值对 BPT 进行剪枝，并输出图像标签
        function label = PruneByHomogeneity(obj, threshold)

            n_leaves = obj.NumLeaves;
            n_nodes = obj.Root;

            % 记录每个结点的左、右孩子ID及其同质度
            us = zeros(n_nodes, 1, "int32");
            vs = zeros(n_nodes, 1, "int32");
            phis = zeros(n_nodes, 1);
            for i = 1:length(obj.Ps)
                p = obj.Ps(i);
                us(p) = obj.Us(i);
                vs(p) = obj.Vs(i);
                phis(p) = obj.Phis(i);
            end

            % 自顶向下剪枝，找到需要剪掉的所有结点ID
            global IPSARTMexHost;
            if isempty(IPSARTMexHost)
                IPSARTMexHost = mexhost();
            end
            pruned = IPSARTMexHost.feval("internal__mex_bridge", ...
                "PolBPT::Prune", us, vs, phis, threshold);

            % 删除生成这些结点的合并记录
            us = obj.Us;
            us(pruned - n_leaves) = [];
            vs = obj.Vs;
            vs(pruned - n_leaves) = [];
            ps = obj.Ps;
            ps(pruned - n_leaves) = [];

            % 剪枝后，利用路径压缩找到每个区域的根结点ID
            map_idx = 1:obj.Root;
            map_idx(us) = ps;
            map_idx(vs) = ps;
            final_idx = IPSARTMexHost.feval("internal__mex_bridge", ...
                "PolBPT::FindRoot", cast(map_idx, 'int32'));

            % 前 1..num_leaves 个结点是每个像素的标签
            pixel_idx = final_idx(1:obj.NumLeaves);

            [~, ~, label] = unique(pixel_idx);
            label = reshape(label, obj.ImageSize);
        end

    end

end
