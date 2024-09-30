function ifgFlt = GoldsteinPhaseFilting(ifg, windowSize, overlap, alphaOrCorr)

arguments
    ifg (:,:)
    windowSize (1,2) {mustBeInteger}
    overlap (1,2) {mustBeInteger}
    alphaOrCorr double
end

if windowSize(1) < 2 * overlap(1)
    error("Overlap in rows must be less than half of window size in rows");
end
row_step = windowSize(1) - overlap(1);

if windowSize(2) < 2 * overlap(2)
    error("Overlap in cols must be less than half of window size in cols");
end
col_step = windowSize(2) - overlap(2);

use_coh_flag = false;
if numel(alphaOrCorr) > 1
    if ~all(size(ifg) == size(alphaOrCorr))
        error("The given correlation coeffcient must have the same size as 'ifg'");
    end
    use_coh_flag = true;
end

[rows,cols] = size(ifg);

n_iter_row = ceil(rows / row_step);

n_iter_col = ceil(cols / col_step);

phase_flt = zeros([rows, cols]);

kernel = fspecial("gaussian", windowSize);

for pat_id_c = 1:n_iter_col
    for pat_id_r = 1:n_iter_row
        r_low = 1 + (pat_id_r - 1) * row_step;
        r_up = min(r_low + windowSize(1) - 1, rows);
        c_low = 1 + (pat_id_c - 1) * col_step;
        c_up = min(c_low + windowSize(2) - 1, cols);

        patch = ifg(r_low:r_up, c_low:c_up);

        z_ifg = fft2(patch);

        smooth = conv2(abs(z_ifg), kernel, "same");
        if use_coh_flag
            alpha = 1 - mean(alphaOrCorr(r_low:r_up, c_low:c_up), "all");
        else
            alpha = alphaOrCorr; 
        end
        smooth = (smooth / max(smooth(:))) .^ alpha;
        phase_flt(r_low:r_up, c_low:c_up) = angle(ifft2(z_ifg .* smooth));
    end
end

ifgFlt = abs(ifg) .* exp(1i * phase_flt);

end