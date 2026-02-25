function ifgFlt = GoldsteinPhaseFilting(ifg, windowSize, overlap, alphaOrCorr)
arguments
    ifg (:,:) {mustBeNumeric}
    windowSize (1,1) {mustBeInteger} = 32 % 建议使用正方形窗口，通常为 32 或 64
    overlap (1,1) {mustBeInteger} = 24    % 建议重叠 75% 以上以消除边界
    alphaOrCorr double = 0.5              % 默认固定 alpha，或传入相干性矩阵
end

[rows, cols] = size(ifg);
% 预计算 Hanning 窗，用于平滑拼合 Patch，消除块状效应
win = hanning(windowSize) * hanning(windowSize)';

% 计算步长
step = windowSize - overlap;
if step <= 0, error("Overlap must be smaller than windowSize"); end

% 准备输出缓冲区（使用复数相加，最后求相位）
acc_ifg = zeros(rows, cols, 'like', ifg);
acc_weight = zeros(rows, cols, 'single');

% 处理相干性矩阵
use_coh = numel(alphaOrCorr) > 1;
if use_coh && ~all(size(alphaOrCorr) == size(ifg))
    error("Size of alphaOrCorr must match ifg");
end

% 预计算格网坐标
r_starts = 1:step:(rows - windowSize + 1);
c_starts = 1:step:(cols - windowSize + 1);

% Goldstein 滤波主循环
for r = r_starts
    r_range = r:(r + windowSize - 1);
    for c = c_starts
        c_range = c:(c + windowSize - 1);

        % 1. 提取 Patch 并变换到频域
        patch = ifg(r_range, c_range);
        patch_fft = fft2(patch);

        % 2. 频域幅度谱平滑 (核心优化：简单的低通平滑)
        mag = abs(patch_fft);
        % 使用简单的均值模糊平滑幅度谱（Goldstein 标准做法）
        mag_smooth = imfilter(mag, ones(3,3)/9, 'circular');

        % 3. 确定自适应参数 alpha
        if use_coh
            % 修正逻辑：coh 高时 alpha 小，coh 低时 alpha 大 (0~1)
            % 这里的均值通常反映局部质量
            avg_coh = mean(alphaOrCorr(r_range, c_range), 'all');
            alpha = 1 - avg_coh;
        else
            alpha = alphaOrCorr;
        end

        % 4. 频域滤波：H(u,v) = S(u,v)^alpha
        % 注意：需要归一化幅度以防止数值爆炸
        mag_op = (mag_smooth ./ (max(mag_smooth(:)) + eps)) .^ alpha;
        filtered_patch_fft = patch_fft .* mag_op;

        % 5. 变换回空域并进行加权堆叠
        filtered_patch = ifft2(filtered_patch_fft);

        % 核心改进：使用复数空间进行 Hanning 加权叠加，消除边缘跳变
        acc_ifg(r_range, c_range) = acc_ifg(r_range, c_range) + filtered_patch .* win;
        acc_weight(r_range, c_range) = acc_weight(r_range, c_range) + win;
    end
end

% 6. 归一化权重并保持原始幅度
% 如果需要保留幅度，使用 original_abs * exp(j * filtered_phase)
% 如果是单纯滤波，直接输出相位
mask = acc_weight > 0;
ifgFlt = zeros(rows, cols, 'like', ifg);
ifgFlt(mask) = abs(ifg(mask)) .* exp(1i * angle(acc_ifg(mask)));

end