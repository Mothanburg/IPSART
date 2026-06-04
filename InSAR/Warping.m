% Warping - 对配准后的图像进行重采样
% 此函数具备一定的可用性
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function result = Warping(slave, gcpMaster, gcpSlave, gcpCorr, polyn_order, method)
arguments
    slave
    gcpMaster (:,2)
    gcpSlave (:,2)
    gcpCorr (:,1) = []
    polyn_order {mustBeMember(polyn_order, [1 2 3 4 5])} = 2
    method string {mustBeMember(method, ["Cubic", "Sinc"])} = "Sinc"
end

% fitting polynomial
[H,W] = size(slave);
[poly_L,poly_P] = fit_polynomial_model( ...
    gcpMaster, gcpSlave, gcpCorr, polyn_order, H, W);

% calculate query
[P_mat,L_mat] = meshgrid(1:W, 1:H);
[Lq,Pq] = apply_polynomial_model( ...
    poly_L, poly_P, L_mat, P_mat, polyn_order, H, W);

% warp
if method == "Sinc"
    result = sinc_interp(slave, Pq, Lq, 8);
elseif method == "Cubic"
    GI = griddedInterpolant({1:H, 1:W}, slave, 'cubic', 'none');
    result = GI(Lq, Pq);
else
    result = [];
end

end


% --- Polynomial fitting ---
function [coeffL, coeffP] = fit_polynomial_model(mst, slv, corr, order, H, W)
% normalize
L_norm = 2 * (mst(:,1) - 1) / (H - 1) - 1;
P_norm = 2 * (mst(:,2) - 1) / (W - 1) - 1;
% fitting
A = construct_vander(L_norm, P_norm, order);
if isempty(corr)
    A_combine = A;
    L_combine = slv(:,1);
    P_combine = slv(:,2);
else
    Q = diag(corr / mean(corr));
    A_combine = A' * Q * A;
    L_combine = A' * Q * slv(:,1);
    P_combine = A' * Q * slv(:,2);
end
coeffL = A_combine \ L_combine;
coeffP = A_combine \ P_combine;
end

function [Lq, Pq] = apply_polynomial_model(cL, cP, L, P, order, H, W)
Ln = 2 * (L(:) - 1) / (H - 1) - 1;
Pn = 2 * (P(:) - 1) / (W - 1) - 1;
A = construct_vander(Ln, Pn, order);
Lq = reshape(A * cL, size(L));
Pq = reshape(A * cP, size(P));
end

function A = construct_vander(L, P, order)
n = length(L);
% generate polynomial terms
num_terms = (order + 1) * (order + 2) / 2;
A = zeros(n, num_terms);
col = 1;
for i = 0:order
    for j = 0:i
        A(:, col) = (L.^(i-j)) .* (P.^j);
        col = col + 1;
    end
end
end


% --- Sinc interpolation ---
function ResultImage = sinc_interp(image, Pq, Lq, kernel_len)
% Calculate Sinc LUT
INTERVAL = 127;
NInterval = INTERVAL + 1;
dx = 1 / INTERVAL;
x_axis = (1 - kernel_len / 2 : kernel_len / 2)';
x_axis_mat = zeros(kernel_len, NInterval);
temp_x = x_axis;
for i = 1 : NInterval
    x_axis_mat(:, i) = temp_x;
    temp_x = temp_x - dx;
end
sinc_mat = sinc(x_axis_mat);
sinc_mat = sinc_mat ./ sum(sinc_mat, 1); % normalize

% 2d interpolation
% horizontal
tmp_res = sinc_resample_1d(image, Pq, sinc_mat, INTERVAL, 'h');
% vertical
ResultImage = sinc_resample_1d(tmp_res, Lq, sinc_mat, INTERVAL, 'v');
end

function out = sinc_resample_1d(img, query, sinc_lut, INTERVAL, mode)
[lines, pixels] = size(img);
out = zeros(lines, pixels, 'like', img);
K = size(sinc_lut, 1); % KernelLength
half_K = K / 2;

% get index
int_idx = floor(query);
frac_idx = query - int_idx;

% find matching kernel
krnl_idx = round(frac_idx * INTERVAL) + 1;
krnl_idx(krnl_idx > (INTERVAL + 1)) = INTERVAL + 1;
krnl_idx(krnl_idx < 1) = 1;

% interp on the whole image
for k = -half_K + 1 : half_K
    if mode == 'h'
        neighbor_idx = min(max(int_idx + k, 1), pixels);
        W = reshape(sinc_lut(k + half_K, krnl_idx(:)), lines, pixels);
        out = out + img((neighbor_idx-1)*lines + (1:lines)') .* W;
    else
        neighbor_idx = min(max(int_idx + k, 1), lines);
        W = reshape(sinc_lut(k + half_K, krnl_idx(:)), lines, pixels);
        out = out + img((0:pixels-1)*lines + neighbor_idx) .* W;
    end
end
end