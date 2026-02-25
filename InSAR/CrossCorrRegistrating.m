function [gcp_mst,gcp_slv,gcp_corr] = CrossCorrRegistrating( ...
    master, slave, gcp_count, opt)
arguments
    master
    slave
    gcp_count (1,1) {mustBeInteger} = 2000
    opt.COARSE_WINDOW_SIZE (1,2) = [128 128]
    opt.FINE_WINDOW_SIZE (1,2) = [32 32]
    opt.FINE_OVERSAMPLE_FACTOR (1,2) = [8 8]
    opt.PNR_THRESHOLD_DB = 10.0
    opt.ENABLE_QUADRATIC_FIT logical = true
    opt.CORR_CALCULAION string {mustBeMember(opt.CORR_CALCULAION, ["Complex" "Magnitude"])} = "Magnitude"
end

if gcp_count <= 0
    gcp_count = ...
        2 * round(numel(master) / prod(window_size));
end

% --- generate ground control points ---
[lines,pixels] = size(master);
num_l = round(sqrt(gcp_count * lines / pixels));
num_p = round(sqrt(gcp_count * pixels / lines));
spacing_l = lines / (num_l + 1);
spacing_p = pixels / (num_p + 1);

ls = (1:num_l) * spacing_l;
ps = (1:num_p) * spacing_p;
[P,L] = meshgrid(ps, ls);

init_gcps = round([L(:), P(:)]);
num_gcp = size(init_gcps, 1);

% --- registration ---
if opt.CORR_CALCULAION == "Complex"
    master_data = parallel.pool.Constant(master);
    slave_data = parallel.pool.Constant(slave);
else
    master_data = parallel.pool.Constant(abs(master));
    slave_data = parallel.pool.Constant(abs(slave));
end

res_slv = zeros(num_gcp, 2);
res_corr = zeros(num_gcp, 1);
parfor idx = 1:num_gcp
    mst_l = init_gcps(idx,1);
    mst_p = init_gcps(idx,2);

    % --- Coarse Coregistration ---
    co_mat = cross_correlation( ...
        master_data, slave_data, ...
        mst_l, mst_p, ...
        mst_l, mst_p, ...
        opt.COARSE_WINDOW_SIZE(1), opt.COARSE_WINDOW_SIZE(2), ...
        [1 1]);
    if isnan(co_mat)
        continue;
    end
    % find pixel-level offset
    [max_val,ind] = max(co_mat, [], "all");
    [cm_h,cm_w] = size(co_mat);
    [r,c] = ind2sub([cm_h, cm_w], ind);
    [dl,dp] = peakpos_to_offset(r, c, cm_h, cm_w);
    % calculate peak-noise ratio (PNR)
    background = co_mat;
    bg_r = max(1, r-1):min(opt.COARSE_WINDOW_SIZE(1), r+1);
    bg_c = max(1, c-1):min(opt.COARSE_WINDOW_SIZE(2), c+1);
    background(bg_r,bg_c) = 0;
    avg_noise = sum(background(:)) / (prod(opt.COARSE_WINDOW_SIZE) - 9);
    pnr_db = 10 * log10(max_val / (avg_noise + eps));
    if pnr_db < opt.PNR_THRESHOLD_DB
        continue;
    end

    % --- Fine Coregistration ---
    slv_l = mst_l + dl;
    slv_p = mst_p + dp;
    co_mat = cross_correlation( ...
        master_data, slave_data, ...
        mst_l, mst_p, ...
        slv_l, slv_p, ...
        opt.FINE_WINDOW_SIZE(1), opt.FINE_WINDOW_SIZE(2), ...
        opt.FINE_OVERSAMPLE_FACTOR);
    if isnan(co_mat)
        continue;
    end
    % find sub-pixel-level offset
    [peak,ind] = max(co_mat, [], "all");
    [cm_h,cm_w] = size(co_mat);
    [r,c] = ind2sub([cm_h, cm_w], ind);
    if opt.ENABLE_QUADRATIC_FIT && r > 1 && r < cm_h && c > 1 && c < cm_w
        % --- quadratic surface fitting ---
        z = co_mat(r-1:r+1, c-1:c+1);
        % quadratic fitting: f(x,y) = ax^2 + by^2 + cxy + dx + ey + f
        % use the least squares method to solve for the 0-derivative point
        % which is the pixel offset
        denom_l = (2*z(2,2) - z(1,2) - z(3,2));
        denom_p = (2*z(2,2) - z(2,1) - z(2,3));
        if denom_l ~= 0 && denom_p ~= 0
            dr = (z(3,2) - z(1,2)) / (2 * denom_l);
            dc = (z(2,3) - z(2,1)) / (2 * denom_p);
            r = r + dr;
            c = c + dc;
        end
    end
    [dl,dp] = peakpos_to_offset(r, c, cm_h, cm_w);
    slv_l = slv_l + dl / opt.FINE_OVERSAMPLE_FACTOR(1);
    slv_p = slv_p + dp / opt.FINE_OVERSAMPLE_FACTOR(2);

    res_slv(idx,:) = [slv_l slv_p];
    res_corr(idx) = peak;
end

% select results
valid = res_corr > 0;
gcp_mst = init_gcps(valid, :);
gcp_slv = res_slv(valid, :);
gcp_corr = res_corr(valid);
end


% --- registration engine ---
function corr_mat = cross_correlation( ...
    master, slave, ...
    mst_l, mst_p, ...
    slv_l, slv_p, ...
    wd_h, wd_w, ...
    ovs_factors)
corr_mat = nan;

% check boundary
[H,W] = size(master.Value);
half_h = wd_h / 2;
half_w = wd_w / 2;
if slv_l - half_h < 1 || slv_l + half_h > H || slv_p - half_w < 1 || slv_p + half_w > W
    return;
end
% master and slave sampling window
m_r = (mst_l - half_h):(mst_l + half_h - 1);
m_c = (mst_p - half_w):(mst_p + half_w - 1);
mst_wd = master.Value(m_r,m_c);
mst_wd = mst_wd - mean(mst_wd(:));
s_r = (slv_l - half_h):(slv_l + half_h - 1);
s_c = (slv_p - half_w):(slv_p + half_w - 1);
slv_wd = slave.Value(s_r,s_c);
slv_wd = slv_wd - mean(slv_wd(:));
% calculate correlation using FFT
corr_spec = fft2(mst_wd) .* conj(fft2(slv_wd));
% upsampling
if any(ovs_factors > 1)
    corr_spec = oversample_spectrum(corr_spec, ovs_factors(1), ovs_factors(2));
end
corr_up = abs(ifft2(corr_spec));
% normalize
corr_mat = corr_up / (norm(mst_wd, "fro") * norm(slv_wd, "fro"));
end

function result = oversample_spectrum(spec, h_factor, w_factor)
[h,w] = size(spec);
h_2 = fix(h / 2);
w_2 = fix(w / 2);
result = zeros(h * h_factor, w * w_factor, 'like', spec);
% spectrum filling
result(1:h_2, 1:w_2) = spec(1:h_2, 1:w_2);
result(1:h_2, (end - (w - w_2) + 1):end) = spec(1:h_2, (w_2 + 1):end);
result((end - (h - h_2) + 1):end, 1:w_2) = spec((h_2 + 1):end, 1:w_2);
result((end - (h - h_2) + 1):end, (end - (w - w_2) + 1):end) = spec((h_2 + 1):end, (w_2 + 1):end);
% keep energy consistent
result = result * (h_factor * w_factor);
end

function [dl,dp] = peakpos_to_offset(pr,pc,h,w)
shift_l = pr - 1;
if shift_l > h / 2
    shift_l = shift_l - h;
end
shift_p = pc - 1;
if shift_p > w / 2
    shift_p = shift_p - w;
end
dl = -shift_l;
dp = -shift_p;
end