% CrossCorrRegistrating - 互相关配准
% 注意！！此函数尚未完善，不建议使用
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [gcp_mst, gcp_slv, gcp_corr] = CrossCorrRegistrating( ...
    master, slave, gcp_count, opt)
arguments
    master
    slave
    gcp_count (1,1) {mustBeInteger} = 2000

    opt.COARSE_WINDOW_SIZE (1,2) double = [128 128]
    opt.FINE_WINDOW_SIZE   (1,2) double = [64 64]

    % V2: fine stage uses local DFT upsampling instead of full-spectrum oversampling
    opt.FINE_UPSAMPLE_FACTOR (1,1) double = 8
    opt.FINE_LOCAL_SPAN_PIX  (1,1) double = 3.0   % local search span in original-pixel units

    opt.COARSE_PSR_THRESHOLD (1,1) double = 5.0
    opt.FINE_PSR_THRESHOLD   (1,1) double = 5.0
    opt.SECOND_PEAK_RATIO_MAX (1,1) double = 0.85

    opt.ENABLE_QUADRATIC_FIT logical = true
    opt.ENABLE_WINDOW logical = true
    opt.ENABLE_ROBUST_OUTLIER_REJECTION logical = true

    % robust outlier filtering params
    opt.OUTLIER_MAD_SCALE (1,1) double = 4.0
    opt.MIN_VALID_GCP_COUNT (1,1) double = 30

    opt.CORR_CALCULAION string ...
        {mustBeMember(opt.CORR_CALCULAION, ["Complex" "Magnitude"])} = "Magnitude"
end

if gcp_count <= 0
    gcp_count = 2 * round(numel(master) / prod(opt.COARSE_WINDOW_SIZE));
end

% --- generate initial GCPs on a grid ---
[lines, pixels] = size(master);
num_l = round(sqrt(gcp_count * lines / pixels));
num_p = round(sqrt(gcp_count * pixels / lines));

spacing_l = lines / (num_l + 1);
spacing_p = pixels / (num_p + 1);

ls = (1:num_l) * spacing_l;
ps = (1:num_p) * spacing_p;
[P, L] = meshgrid(ps, ls);

init_gcps = round([L(:), P(:)]);
num_gcp = size(init_gcps, 1);

% --- choose data mode ---
if opt.CORR_CALCULAION == "Complex"
    master_data = parallel.pool.Constant(master);
    slave_data  = parallel.pool.Constant(slave);
else
    master_data = parallel.pool.Constant(abs(master));
    slave_data  = parallel.pool.Constant(abs(slave));
end

res_slv   = nan(num_gcp, 2);
res_corr  = nan(num_gcp, 1);
res_psr   = nan(num_gcp, 1);
res_shift = nan(num_gcp, 2);

parfor idx = 1:num_gcp
    mst_l = init_gcps(idx,1);
    mst_p = init_gcps(idx,2);

    % =========================================================
    % 1) coarse registration: linear correlation
    % =========================================================
    [co_mat, ~] = cross_correlation_linear( ...
        master_data, slave_data, ...
        mst_l, mst_p, ...
        mst_l, mst_p, ...
        opt.COARSE_WINDOW_SIZE(1), opt.COARSE_WINDOW_SIZE(2), ...
        opt.ENABLE_WINDOW);

    if isscalar(co_mat) && isnan(co_mat)
        continue;
    end

    [max_val, ind] = max(co_mat, [], "all");
    [cm_h, cm_w] = size(co_mat);
    [r, c] = ind2sub([cm_h, cm_w], ind);

    q1 = peak_quality(co_mat, [r, c]);
    if ~q1.valid
        continue;
    end
    if q1.psr < opt.COARSE_PSR_THRESHOLD
        continue;
    end
    if q1.second_ratio > opt.SECOND_PEAK_RATIO_MAX
        continue;
    end

    [dl0, dp0] = peakpos_to_offset(r, c, cm_h, cm_w);

    slv_l0 = mst_l + dl0;
    slv_p0 = mst_p + dp0;

    % =========================================================
    % 2) fine registration (V2):
    %    linear correlation -> integer peak -> local DFT upsampling
    % =========================================================
    [fine_mat, fine_spec] = cross_correlation_linear( ...
        master_data, slave_data, ...
        mst_l, mst_p, ...
        slv_l0, slv_p0, ...
        opt.FINE_WINDOW_SIZE(1), opt.FINE_WINDOW_SIZE(2), ...
        opt.ENABLE_WINDOW);

    if isscalar(fine_mat) && isnan(fine_mat)
        continue;
    end

    [peak0, ind0] = max(fine_mat, [], "all");
    [fh, fw] = size(fine_mat);
    [r0, c0] = ind2sub([fh, fw], ind0);

    q2 = peak_quality(fine_mat, [r0, c0]);
    if ~q2.valid
        continue;
    end
    if q2.psr < opt.FINE_PSR_THRESHOLD
        continue;
    end
    if q2.second_ratio > opt.SECOND_PEAK_RATIO_MAX
        continue;
    end

    [dl_int, dp_int] = peakpos_to_offset(r0, c0, fh, fw);

    % local DFT upsampling around integer peak
    Cup = local_dft_upsample( ...
        fine_spec, ...
        [dl_int, dp_int], ...
        opt.FINE_UPSAMPLE_FACTOR, ...
        opt.FINE_LOCAL_SPAN_PIX);

    Sup = abs(Cup);
    [peak, indu] = max(Sup, [], "all");
    [uh, uw] = size(Sup);
    [ru, cu] = ind2sub([uh, uw], indu);

    q3 = peak_quality(Sup, [ru, cu]);
    if ~q3.valid
        continue;
    end
    if q3.psr < opt.FINE_PSR_THRESHOLD
        continue;
    end
    if q3.second_ratio > opt.SECOND_PEAK_RATIO_MAX
        continue;
    end

    [dlu, dpu] = local_peakpos_to_offset( ...
        ru, cu, [uh, uw], ...
        opt.FINE_UPSAMPLE_FACTOR);

    % optional final 2D quadratic tweak on upsampled local surface
    if opt.ENABLE_QUADRATIC_FIT && has_3x3([uh, uw], [ru, cu])
        z = Sup(ru-1:ru+1, cu-1:cu+1);
        [edl, edp, ok] = fit_quadratic_2d(z);
        if ok
            dlu = dlu + edl / opt.FINE_UPSAMPLE_FACTOR;
            dpu = dpu + edp / opt.FINE_UPSAMPLE_FACTOR;
        end
    end

    slv_l = slv_l0 + dl_int + dlu;
    slv_p = slv_p0 + dp_int + dpu;

    res_slv(idx,:)   = [slv_l, slv_p];
    res_corr(idx,1)  = peak;
    res_psr(idx,1)   = q3.psr;
    res_shift(idx,:) = [slv_l - mst_l, slv_p - mst_p];
end

% ============================================================
% 3) collect valid results
% ============================================================
valid = all(isfinite(res_slv), 2) & isfinite(res_corr) & isfinite(res_psr);
gcp_mst  = init_gcps(valid, :);
gcp_slv  = res_slv(valid, :);
gcp_corr = res_corr(valid);
gcp_psr  = res_psr(valid);
gcp_shift = res_shift(valid, :);

% ============================================================
% 4) robust outlier rejection
% ============================================================
if opt.ENABLE_ROBUST_OUTLIER_REJECTION && size(gcp_shift,1) >= opt.MIN_VALID_GCP_COUNT
    keep = robust_shift_filter(gcp_shift, gcp_psr, opt.OUTLIER_MAD_SCALE);

    gcp_mst  = gcp_mst(keep, :);
    gcp_slv  = gcp_slv(keep, :);
    gcp_corr = gcp_corr(keep, :);
end

end


% ============================================================
% LINEAR cross-correlation
% returns:
%   corr_mat : abs(fftshift(ifft2(FA .* conj(FB))))
%   corr_spec: FA .* conj(FB), on linear-correlation grid size
% ============================================================
function [corr_mat, corr_spec] = cross_correlation_linear( ...
    master, slave, ...
    mst_l, mst_p, ...
    slv_l, slv_p, ...
    wd_h, wd_w, ...
    enable_window)

corr_mat  = nan;
corr_spec = nan;

[H, W] = size(master.Value);
half_h = wd_h / 2;
half_w = wd_w / 2;

% boundary check: BOTH master and slave
if mst_l - half_h < 1 || mst_l + half_h > H || ...
        mst_p - half_w < 1 || mst_p + half_w > W
    return;
end
if slv_l - half_h < 1 || slv_l + half_h > H || ...
        slv_p - half_w < 1 || slv_p + half_w > W
    return;
end

% sample windows
m_r = (mst_l - half_h):(mst_l + half_h - 1);
m_c = (mst_p - half_w):(mst_p + half_w - 1);
s_r = (slv_l - half_h):(slv_l + half_h - 1);
s_c = (slv_p - half_w):(slv_p + half_w - 1);

mst_wd = master.Value(m_r, m_c);
slv_wd = slave.Value(s_r, s_c);

% remove DC
mst_wd = mst_wd - mean(mst_wd(:));
slv_wd = slv_wd - mean(slv_wd(:));

% window
if enable_window
    wr = hann(wd_h);
    wc = hann(wd_w);
    W2 = wr * wc.';
    mst_wd = mst_wd .* W2;
    slv_wd = slv_wd .* W2;
end

nm = norm(mst_wd, "fro");
ns = norm(slv_wd, "fro");
if nm < eps || ns < eps
    return;
end

P = 2 * wd_h - 1;
Q = 2 * wd_w - 1;

FA = fft2(mst_wd, P, Q);
FB = fft2(slv_wd, P, Q);

corr_spec = FA .* conj(FB);

corr = fftshift(ifft2(corr_spec));
corr_mat = abs(corr) / (nm * ns + eps);
end


% ============================================================
% local DFT upsampling around coarse integer peak
%
% corr_spec : frequency product on linear-correlation grid
% coarse_shift = [dl_int, dp_int], in original pixels
% up         : scalar upsampling factor
% spanPix    : local window size in original-pixel units
%
% returns a local refined correlation surface only around the peak
% ============================================================
function Cup = local_dft_upsample(corr_spec, coarse_shift, up, spanPix)

[P, Q] = size(corr_spec);

halfSpan = spanPix / 2;
u = coarse_shift(1) + (-halfSpan : 1/up : halfSpan);   % line shift grid
v = coarse_shift(2) + (-halfSpan : 1/up : halfSpan);   % pixel shift grid

% frequency indices consistent with fftshifted spectrum
m = ifftshift((-floor(P/2):ceil(P/2)-1)).';
n = ifftshift((-floor(Q/2):ceil(Q/2)-1));

S = fftshift(corr_spec);

Er = exp(1i * 2*pi/P * (m * u));      % P x Nu
Ec = exp(1i * 2*pi/Q * (v(:) * n));   % Nv x Q

Cup = (Er' * S * Ec.') / (P * Q);     % Nu x Nv
end


% ============================================================
% convert local upsampled peak index to subpixel offset
% relative to the local integer-peak center
% ============================================================
function [dl, dp] = local_peakpos_to_offset(pr, pc, sz, up)
center_r = (sz(1) + 1) / 2;
center_c = (sz(2) + 1) / 2;

dl = (pr - center_r) / up;
dp = (pc - center_c) / up;
end


% ============================================================
% peak index to centered offset for a fftshift()'ed linear correlation map
% ============================================================
function [dl, dp] = peakpos_to_offset(pr, pc, h, w)
dl = pr - ceil(h / 2);
dp = pc - ceil(w / 2);
end


% ============================================================
% peak quality using PSR + second peak ratio
% ============================================================
function q = peak_quality(corr_mat, idx)

r = idx(1);
c = idx(2);
[h, w] = size(corr_mat);

if r < 1 || r > h || c < 1 || c > w
    q.valid = false;
    q.psr = -inf;
    q.second_ratio = inf;
    q.main_peak = NaN;
    q.second_peak = NaN;
    return;
end

main_peak = corr_mat(r, c);

mask = true(h, w);
rr = max(1, r-1):min(h, r+1);
cc = max(1, c-1):min(w, c+1);
mask(rr, cc) = false;

side = corr_mat(mask);
if isempty(side)
    q.valid = false;
    q.psr = -inf;
    q.second_ratio = inf;
    q.main_peak = main_peak;
    q.second_peak = NaN;
    return;
end

side_mean = mean(side);
side_std  = std(side);
second_peak = max(side);

psr = (main_peak - side_mean) / (side_std + eps);
second_ratio = second_peak / (main_peak + eps);

q.valid = isfinite(psr) && isfinite(second_ratio) && isfinite(main_peak);
q.psr = psr;
q.second_ratio = second_ratio;
q.main_peak = main_peak;
q.second_peak = second_peak;
end


% ============================================================
% true 2D quadratic fitting on 3x3 patch
% z(x,y)=a*x^2+b*y^2+c*x*y+d*x+e*y+f
% ============================================================
function [dr, dc, ok] = fit_quadratic_2d(z)

[RR, CC] = ndgrid(-1:1, -1:1);
A = [RR(:).^2, CC(:).^2, RR(:).*CC(:), RR(:), CC(:), ones(9,1)];
coef = A \ z(:);

a = coef(1);
b = coef(2);
c = coef(3);
d = coef(4);
e = coef(5);

H = [2*a, c;
    c, 2*b];
g = [d; e];

ok = rcond(H) > 1e-10 && all(eig(H) < 0);
if ~ok
    dr = 0;
    dc = 0;
    return;
end

v = -H \ g;
dr = v(1);
dc = v(2);

if max(abs(v)) > 1.0
    ok = false;
    dr = 0;
    dc = 0;
end
end


% ============================================================
% robust outlier filtering based on global shift median + MAD
% PSR is used as a soft preference when MAD becomes tiny
% ============================================================
function keep = robust_shift_filter(shift, psr, mad_scale)

dl = shift(:,1);
dp = shift(:,2);

med_dl = median(dl, "omitnan");
med_dp = median(dp, "omitnan");

mad_dl = median(abs(dl - med_dl), "omitnan");
mad_dp = median(abs(dp - med_dp), "omitnan");

% convert MAD to robust sigma estimate
sig_dl = 1.4826 * max(mad_dl, eps);
sig_dp = 1.4826 * max(mad_dp, eps);

z_dl = abs(dl - med_dl) / sig_dl;
z_dp = abs(dp - med_dp) / sig_dp;

keep = (z_dl <= mad_scale) & (z_dp <= mad_scale);

% if too strict due to very concentrated offsets, relax by PSR ranking
if nnz(keep) < max(20, round(0.2 * numel(keep)))
    score = psr;
    score(~isfinite(score)) = -inf;

    [~, ord] = sort(score, "descend");
    keep = false(size(keep));
    keep(ord(1:min(numel(ord), max(20, round(0.5 * numel(ord)))))) = true;
end
end


% ============================================================
function tf = has_3x3(sz, idx)
h = sz(1); w = sz(2);
r = idx(1); c = idx(2);
tf = (r > 1) && (r < h) && (c > 1) && (c < w);
end

% function [gcp_mst,gcp_slv,gcp_corr] = CrossCorrRegistrating( ...
%     master, slave, gcp_count, opt)
% arguments
%     master
%     slave
%     gcp_count (1,1) {mustBeInteger} = 2000
%     opt.COARSE_WINDOW_SIZE (1,2) = [128 128]
%     opt.FINE_WINDOW_SIZE (1,2) = [32 32]
%     opt.FINE_OVERSAMPLE_FACTOR (1,2) = [8 8]
%     opt.PNR_THRESHOLD_DB = 10.0
%     opt.ENABLE_QUADRATIC_FIT logical = true
%     opt.CORR_CALCULAION string {mustBeMember(opt.CORR_CALCULAION, ["Complex" "Magnitude"])} = "Magnitude"
% end
%
% if gcp_count <= 0
%     gcp_count = ...
%         2 * round(numel(master) / prod(window_size));
% end
%
% % --- generate ground control points ---
% [lines,pixels] = size(master);
% num_l = round(sqrt(gcp_count * lines / pixels));
% num_p = round(sqrt(gcp_count * pixels / lines));
% spacing_l = lines / (num_l + 1);
% spacing_p = pixels / (num_p + 1);
%
% ls = (1:num_l) * spacing_l;
% ps = (1:num_p) * spacing_p;
% [P,L] = meshgrid(ps, ls);
%
% init_gcps = round([L(:), P(:)]);
% num_gcp = size(init_gcps, 1);
%
% % --- registration ---
% if opt.CORR_CALCULAION == "Complex"
%     master_data = parallel.pool.Constant(master);
%     slave_data = parallel.pool.Constant(slave);
% else
%     master_data = parallel.pool.Constant(abs(master));
%     slave_data = parallel.pool.Constant(abs(slave));
% end
%
% res_slv = zeros(num_gcp, 2);
% res_corr = zeros(num_gcp, 1);
% parfor idx = 1:num_gcp
%     mst_l = init_gcps(idx,1);
%     mst_p = init_gcps(idx,2);
%
%     % --- Coarse Coregistration ---
%     co_mat = cross_correlation( ...
%         master_data, slave_data, ...
%         mst_l, mst_p, ...
%         mst_l, mst_p, ...
%         opt.COARSE_WINDOW_SIZE(1), opt.COARSE_WINDOW_SIZE(2), ...
%         [1 1]);
%     if isnan(co_mat)
%         continue;
%     end
%     % find pixel-level offset
%     [max_val,ind] = max(co_mat, [], "all");
%     [cm_h,cm_w] = size(co_mat);
%     [r,c] = ind2sub([cm_h, cm_w], ind);
%     [dl,dp] = peakpos_to_offset(r, c, cm_h, cm_w);
%     % calculate peak-noise ratio (PNR)
%     background = co_mat;
%     bg_r = max(1, r-1):min(opt.COARSE_WINDOW_SIZE(1), r+1);
%     bg_c = max(1, c-1):min(opt.COARSE_WINDOW_SIZE(2), c+1);
%     background(bg_r,bg_c) = 0;
%     avg_noise = sum(background(:)) / (prod(opt.COARSE_WINDOW_SIZE) - 9);
%     pnr_db = 10 * log10(max_val / (avg_noise + eps));
%     if pnr_db < opt.PNR_THRESHOLD_DB
%         continue;
%     end
%
%     % --- Fine Coregistration ---
%     slv_l = mst_l + dl;
%     slv_p = mst_p + dp;
%     co_mat = cross_correlation( ...
%         master_data, slave_data, ...
%         mst_l, mst_p, ...
%         slv_l, slv_p, ...
%         opt.FINE_WINDOW_SIZE(1), opt.FINE_WINDOW_SIZE(2), ...
%         opt.FINE_OVERSAMPLE_FACTOR);
%     if isnan(co_mat)
%         continue;
%     end
%     % find sub-pixel-level offset
%     [peak,ind] = max(co_mat, [], "all");
%     [cm_h,cm_w] = size(co_mat);
%     [r,c] = ind2sub([cm_h, cm_w], ind);
%     if opt.ENABLE_QUADRATIC_FIT && r > 1 && r < cm_h && c > 1 && c < cm_w
%         % --- quadratic surface fitting ---
%         z = co_mat(r-1:r+1, c-1:c+1);
%         % quadratic fitting: f(x,y) = ax^2 + by^2 + cxy + dx + ey + f
%         % use the least squares method to solve for the 0-derivative point
%         % which is the pixel offset
%         denom_l = (2*z(2,2) - z(1,2) - z(3,2));
%         denom_p = (2*z(2,2) - z(2,1) - z(2,3));
%         if denom_l ~= 0 && denom_p ~= 0
%             dr = (z(3,2) - z(1,2)) / (2 * denom_l);
%             dc = (z(2,3) - z(2,1)) / (2 * denom_p);
%             r = r + dr;
%             c = c + dc;
%         end
%     end
%     [dl,dp] = peakpos_to_offset(r, c, cm_h, cm_w);
%     slv_l = slv_l + dl / opt.FINE_OVERSAMPLE_FACTOR(1);
%     slv_p = slv_p + dp / opt.FINE_OVERSAMPLE_FACTOR(2);
%
%     res_slv(idx,:) = [slv_l slv_p];
%     res_corr(idx) = peak;
% end
%
% % select results
% valid = res_corr > 0;
% gcp_mst = init_gcps(valid, :);
% gcp_slv = res_slv(valid, :);
% gcp_corr = res_corr(valid);
% end
%
%
% % --- registration engine ---
% function corr_mat = cross_correlation( ...
%     master, slave, ...
%     mst_l, mst_p, ...
%     slv_l, slv_p, ...
%     wd_h, wd_w, ...
%     ovs_factors)
% corr_mat = nan;
%
% % check boundary
% [H,W] = size(master.Value);
% half_h = wd_h / 2;
% half_w = wd_w / 2;
% if slv_l - half_h < 1 || slv_l + half_h > H || slv_p - half_w < 1 || slv_p + half_w > W
%     return;
% end
% % master and slave sampling window
% m_r = (mst_l - half_h):(mst_l + half_h - 1);
% m_c = (mst_p - half_w):(mst_p + half_w - 1);
% mst_wd = master.Value(m_r,m_c);
% mst_wd = mst_wd - mean(mst_wd(:));
% s_r = (slv_l - half_h):(slv_l + half_h - 1);
% s_c = (slv_p - half_w):(slv_p + half_w - 1);
% slv_wd = slave.Value(s_r,s_c);
% slv_wd = slv_wd - mean(slv_wd(:));
% % calculate correlation using FFT
% corr_spec = fft2(mst_wd) .* conj(fft2(slv_wd));
% % upsampling
% if any(ovs_factors > 1)
%     corr_spec = oversample_spectrum(corr_spec, ovs_factors(1), ovs_factors(2));
% end
% corr_up = abs(ifft2(corr_spec));
% % normalize
% corr_mat = corr_up / (norm(mst_wd, "fro") * norm(slv_wd, "fro"));
% end
%
% function result = oversample_spectrum(spec, h_factor, w_factor)
% [h,w] = size(spec);
% h_2 = fix(h / 2);
% w_2 = fix(w / 2);
% result = zeros(h * h_factor, w * w_factor, 'like', spec);
% % spectrum filling
% result(1:h_2, 1:w_2) = spec(1:h_2, 1:w_2);
% result(1:h_2, (end - (w - w_2) + 1):end) = spec(1:h_2, (w_2 + 1):end);
% result((end - (h - h_2) + 1):end, 1:w_2) = spec((h_2 + 1):end, 1:w_2);
% result((end - (h - h_2) + 1):end, (end - (w - w_2) + 1):end) = spec((h_2 + 1):end, (w_2 + 1):end);
% % keep energy consistent
% result = result * (h_factor * w_factor);
% end
%
% function [dl,dp] = peakpos_to_offset(pr,pc,h,w)
% shift_l = pr - 1;
% if shift_l > h / 2
%     shift_l = shift_l - h;
% end
% shift_p = pc - 1;
% if shift_p > w / 2
%     shift_p = shift_p - w;
% end
% dl = -shift_l;
% dp = -shift_p;
% end