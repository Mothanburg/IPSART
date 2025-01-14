function result = Warping(master, slave, gcpMaster, gcpSlave, gcpCorr, options)

arguments
    master
    slave
    gcpMaster
    gcpSlave
    gcpCorr = []
    options.WARP_POLYNOMIAL_ORDER = 3
    options.WARP_SOLVING_METHOD = "LS"
    options.RESAMPLE_ESTIMATE_DOPPLER_CENTROID = false
    options.RESAMPLE_INTERPOLATION_METHOD = "16-point sinc"
end

[lines,pixels] = size(slave);
[n_gcp,~] = size(gcpMaster);


% normalize master coordinates into [-2 2] for stability
order = options.WARP_POLYNOMIAL_ORDER;

mst_lines_norm = normalize2(gcpMaster(:,1), 1, lines);
mst_pixels_norm = normalize2(gcpMaster(:,2), 1, pixels);
A = create_polyn_param_mat( ...
    mst_lines_norm, ...
    mst_pixels_norm, ...
    order ...
    );

if ~isempty(gcpCorr)
    Q = diag(gcpCorr / mean(gcpCorr)); % weighting matrix
else
    Q = eye(n_gcp);
end

switch (options.WARP_SOLVING_METHOD)
    case "LS"
        warp_coeff_l = (A' * Q * A) \ ...
            (A' * Q * (gcpSlave(:,1) - gcpMaster(:,1)));
        warp_coeff_p = (A' * Q * A) \ ...
            (A' * Q * (gcpSlave(:,2) - gcpMaster(:,2)));
    otherwise
        error("Unknown error")
end


% Get quering coordinates in slave image
[P,L] = meshgrid(1:pixels, 1:lines);

norm_lines = normalize2(L(:), 1, lines);
norm_pixels = normalize2(P(:), 1, pixels);

dL = apply_polyn(warp_coeff_l, norm_lines, norm_pixels, order);
dL = reshape(dL, size(L));

dP = apply_polyn(warp_coeff_p, norm_lines, norm_pixels, order);
dP = reshape(dP, size(P));

Pq = P + dP;
Lq = L + dL;


% Resampling
interp_method = lower(options.RESAMPLE_INTERPOLATION_METHOD);
switch interp_method
    case {"nearest", "linear"}
        result = interp2( ...
            P, L, ...
            slave, ...
            Pq, Lq, ...
            interp_method, ...
            nan ...
            );
    case "6-point cubic"
        krnl_lut = construct_interp_table(@cc_6p, 6);
        result = interp_image(slave, Lq, Pq, krnl_lut, nan);
    case "8-point sinc"
        krnl_lut = construct_interp_table(@sinc_8p, 8);
        result = interp_image(slave, Lq, Pq, krnl_lut, nan);
    case "16-point sinc"
        krnl_lut = construct_interp_table(@sinc_16p, 16);
        result = interp_image(slave, Lq, Pq, krnl_lut, nan);
    otherwise
        error("Unknown error");
end


end



%----------------- Normalize coordinates into [-2 2] -----------------%
function result = normalize2(value, low, high)
    result = 4 * (value - low) ./ (high - low) - 2;
end



%----------------- Create Polynomial Parameter Matrix -----------------%
function A = create_polyn_param_mat(lines, pixels, order)

len = length(lines);
A = zeros(len, ((order + 1)^2 + order + 1) / 2);
idx = 1;
for i = 0:order
    for j = 0:i
        A(:,idx) = lines.^(i - j) .* pixels.^(j);
        idx = idx + 1;
    end
end

end

function results = apply_polyn(coeff, lines, pixels, order)

A = create_polyn_param_mat(lines, pixels, order);
results = A * coeff;

end



%----------------- Construct Interpolation LUT -----------------%
function [value,axis] = construct_interp_table(krnlFunc, krnlLen)

n_quant = 128;
spacing = 1 / (n_quant - 1);
base_axis = (1:krnlLen) - krnlLen / 2;

value = zeros(n_quant, krnlLen);
axis = zeros(n_quant, krnlLen);
for idx = 1:n_quant
    krnl_axis = base_axis - spacing * (idx - 1);
    axis(idx,:) = krnl_axis;
    value(idx,:) = krnlFunc(krnl_axis);
end

end

function value = sinc_8p(axis)

value = sinc(axis);
value = value / sum(value);

end

function value = sinc_16p(axis)

value = sinc(axis);
value = value / sum(value);

end

function value = cc_6p(axis)

alpha = -0.5;
beta = 0.5;
value = zeros(1, 6);
for i = 1:6
    x = abs(axis(i));
    if x < 1
        value(i) = (alpha - beta + 2) * x^3 - (alpha - beta + 3) * x^2 + 1;
    elseif x < 2
        value(i) = alpha * x^3 - (5 * alpha - beta) * x^2 + ...
            (8 * alpha - 3 * beta) * x^2 - (4 * alpha - 2 * beta);
    elseif x < 3
        value(i) = beta * x^3 - 8 * beta * x^2 + 21 * beta * x^2 - 18 * beta;
    else
        value(i) = 0;
    end
end
value = value / sum(value);

end



%----------------- Image Interpolation -----------------%
function result = interp_image(image, Lq, Pq, lut, extraVal)

result = zeros(size(Lq));
[lines,pixels] = size(Lq);
len = lines * pixels;

[n_quant,krnl_len] = size(lut);

margin = round(krnl_len / 2);
margined = zeros(lines + krnl_len, pixels + krnl_len);
margined(margin + (1:lines),margin + (1:pixels)) = image;

Lq = parallel.pool.Constant(Lq);
Pq = parallel.pool.Constant(Pq);
margined = parallel.pool.Constant(margined);
parfor idx = 1:len
    lq = Lq.Value(idx);
    pq = Pq.Value(idx);

    if lq < 1 || lq > lines || pq < 1 || pq > pixels
        result(idx) = extraVal;
        continue;
    end

    lq_int = floor(lq);
    lq_dec = lq - lq_int;
    krnl_idx = round(lq_dec * (n_quant - 1)) + 1; % which kernel is we want in LUT
    interp_ls = (-krnl_len / 2 + 1 : krnl_len / 2) + margin + lq_int;
    krnl_l = lut(krnl_idx,:);

    pq_int = floor(pq);
    pq_dec = pq - pq_int;
    krnl_idx = round(pq_dec * (n_quant - 1)) + 1;
    interp_ps = (-krnl_len / 2 + 1 : krnl_len / 2) + margin + pq_int;
    krnl_p = lut(krnl_idx,:);

    wd = margined.Value(interp_ls,interp_ps);
    result(idx) = krnl_l * wd * krnl_p';

end


end













