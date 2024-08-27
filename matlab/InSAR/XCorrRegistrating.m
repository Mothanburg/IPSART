function [gcpMst,gcpSlv,gcpCorr] = XCorrRegistrating( ...
    master, ...
    slave, ...
    options ...
    )

arguments
    master
    slave
    options.GCP_NUM = 2000
    % coarse registration
    options.COARSE_WINDOW_SIZE = [128 128]
    options.COARSE_INTERP_FACTOR = [4 4]
    options.COARSE_MAX_ITER = 10
    options.COARSE_GCP_TOLERANCE = 0.25
    options.ENABLE_FINE_REGISTRATION = true

end

% point (l,p) in master + slaveOffset == point (l2,p2) in slave
assert(all(size(master) == size(slave)), ...
    "The image pair must have the same size.");

% Coarse registration
[gcpMst,gcpSlv,gcpCorr] = coarse_registrating( ...
    abs(master), abs(slave), ...
    options.GCP_NUM, ...
    options.COARSE_WINDOW_SIZE(1), options.COARSE_WINDOW_SIZE(2), ...
    options.COARSE_INTERP_FACTOR(1), options.COARSE_INTERP_FACTOR(2), ...
    options.COARSE_GCP_TOLERANCE, ...
    options.COARSE_MAX_ITER ...
    );

% Fine registration
if options.ENABLE_FINE_REGISTRATION

end

end


%%%%%%%%%%%%%%%% Generate GCP %%%%%%%%%%%%%%%% 
function gcps = gen_gcp(lines, pixels, nGCP)

if nGCP >= 0.5 * lines * pixels
    error("The number of GCPs must be greatly less than number of pixels");
end

num_l = round(sqrt(nGCP * lines / pixels));
num_p = round(sqrt(nGCP * pixels / lines));
total_len = num_p * num_l;

spacing_l = lines / (num_l + 1);
spacing_p = pixels / (num_p + 1);

gcps = zeros([total_len 2]);

for ip = 1:num_p
    for il = 1:num_l
        l = round(il * spacing_l);
        p = round(ip * spacing_p);
        idx = sub2ind([num_l, num_p], il, ip);
        gcps(idx,:) = [l, p];
    end
end

end


%%%%%%%%%%%%%%%% Upsampling spectrum by DFT %%%%%%%%%%%%%%%%
function result = spec_upsampling(spec, hFactor, wFactor)

[h,w] = size(spec);
h_2 = h / 2;
w_2 = w / 2;

result = zeros(h * hFactor, w * wFactor);
result(1:h_2,1:w_2) = spec(1:h_2,1:w_2);
result(1:h_2,(end - w_2 + 1):end) = spec(1:h_2,(w_2 + 1):end);
result((end - h_2 + 1):end,1:w_2) = spec((h_2 + 1):end,1:w_2);
result((end - h_2 + 1):end,(end - w_2 + 1):end) = ...
    spec((h_2 + 1):end,(w_2 + 1):end);
% result = result * hFactor * wFactor;

end


%%%%%%%%%%%%%%%% Coarse Registration to generate GCPs %%%%%%%%%%%%%%%%
function [gcp_mst,gcp_slv,gcp_corr] = coarse_registrating(...
    master, slave, ...
    num_gcp, cors_wd_h ,cors_wd_w, ...
    cors_intp_factor_l, cors_intp_factor_p , ...
    gcp_tol, max_iter)

[lines,pixels] = size(master);

init_gcps = gen_gcp(lines, pixels, num_gcp);
[num_gcp,~] = size(init_gcps);

gcp_slv = zeros(size(init_gcps));
gcp_corr = zeros(numel(gcp_slv) / 2);

init_gcps = parallel.pool.Constant(init_gcps);
master = parallel.pool.Constant(master);
slave = parallel.pool.Constant(slave);

parfor idx = 1:num_gcp
    l0 = init_gcps.Value(idx,1);
    p0 = init_gcps.Value(idx,2);

    l = l0;
    p = p0;
    gcp_shift = zeros(1, 2);
    peak = 0;
    for tmp_ = 1:max_iter

        try
            mst_wd = master.Value(...
                (l0 - cors_wd_h / 2):(l0 + cors_wd_h / 2 - 1), ...
                (p0 - cors_wd_w / 2):(p0 + cors_wd_w / 2 - 1) ...
                );
            % because the offset may be decimal, we need interpolation
            slv_wd = zeros(cors_wd_h, cors_wd_w);
            l_start = l - cors_wd_h / 2;
            p_start = p - cors_wd_w / 2;
            for ip = 1:cors_wd_w 
                for il = 1:cors_wd_h
                    cur_p = p_start + ip - 1;
                    cur_l = l_start + il - 1;
                    x00 = slave.Value(floor(cur_l), floor(cur_p));
                    x01 = slave.Value(floor(cur_l), ceil(cur_p));
                    x10 = slave.Value(ceil(cur_l), floor(cur_p));
                    x11 = slave.Value(ceil(cur_l), ceil(cur_p));
                    dp = cur_p - floor(cur_p);
                    dl = cur_l - floor(cur_l);
                    slv_wd(il,ip) = ...
                        x00 + ...
                        dl * (x10 - x00) + ...
                        dp * (x01 - x00) + ...
                        dl * dp * (x11 + x00 - x01 - x10);
                end
            end

            if any(isinf(mst_wd(:))) || any(isnan(mst_wd(:))) || ...
               any(isinf(slv_wd(:))) || any(isnan(slv_wd(:)))
                error("invalid patch");
            end
        catch % encounter accessing error
            l = -1;
            p = -1;
            break; % drop this GCP
        end

        % Calculate XCorr
        mst_wd = mst_wd - mean(mst_wd(:));
        slv_wd = slv_wd - mean(slv_wd(:));

        cor_spec = fft2(mst_wd) .* conj(fft2(slv_wd));
        cor_spec_up = spec_upsampling(...
            cor_spec, ...
            cors_intp_factor_l, ...
            cors_intp_factor_p ...
            );
        cor_up = abs(ifft2(cor_spec_up));

        [peak,ind] = max(cor_up(:));
        [up_h,up_w] = size(cor_up);
        [peak_r,peak_c] = ind2sub([up_h up_w], ind);
        peak_r = peak_r - 1;
        peak_c = peak_c - 1;

        if peak_r > up_h / 2
            gcp_shift(1) = (up_h - peak_r) / cors_intp_factor_l;
        else
            gcp_shift(1) = -peak_r / cors_intp_factor_l;
        end

        if peak_c > up_w / 2
            gcp_shift(2) = (up_w - peak_c) / cors_intp_factor_p;
        else
            gcp_shift(2) = -peak_c / cors_intp_factor_p;
        end

        l = l + gcp_shift(1);
        p = p + gcp_shift(2);

        if all(abs(gcp_shift) <= gcp_tol)
            break;
        end

    end
    
    % drop unstable gcps
    if any(abs(gcp_shift) > gcp_tol)
        l = -1;
        p = -1;
    end

    gcp_slv(idx,:) = [l p];
    gcp_corr(idx) = peak;

end

indices = gcp_slv(:,1) > 0;

gcp_mst = init_gcps.Value(indices,:);
gcp_slv = gcp_slv(indices,:);
gcp_corr = normalize(gcp_corr(indices), "range");

end









