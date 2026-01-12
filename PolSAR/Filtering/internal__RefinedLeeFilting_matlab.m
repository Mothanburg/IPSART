function result = internal__RefinedLeeFilting_matlab(M, lookNum)

height = M.Height;
width = M.Width;
span = M.SPAN;

% Gradient templates
w1 = [-1,0,1;-1,0,1;-1,0,1];
w2 = [0,1,1;-1,0,1;-1,-1,0];
w3 = [1,1,1;0,0,0;-1,-1,-1];
w4 = [1,1,0;1,0,-1;0,-1,-1];
W = cat(3, w1, w2, w3, w4, fliplr(w1), w2', flipud(w3), flipud(w2));

% Using the dilated convolution to get the best Prewitt mask
grads = zeros(height, width, 8);
avg_krnl = ones(3) / 9;
for k = 1:8
    w_k = zeros(5);
    w_k(1:2:5, 1:2:5) = W(:,:,k);

    conv_krnl = conv2(w_k, avg_krnl, 'full');

    grads = imfilter(span, conv_krnl, "replicate", "same");
end

[~, pw_id] = max(grads, [], 3);

% Prewitt masks
pw1 = repmat([0,0,0,1,1,1,1], [7,1]);
pw2 = [1,1,1,1,1,1,1;
    0,1,1,1,1,1,1;
    0,0,1,1,1,1,1;
    0,0,0,1,1,1,1;
    0,0,0,0,1,1,1;
    0,0,0,0,0,1,1;
    0,0,0,0,0,0,1];
pw3 = repmat([1,1,1,1,0,0,0]', [1,7]);
pw4 = flip(pw2, 2);
pw5 = flip(pw1, 2);
pw6 = flip(pw4, 1);
pw7 = flip(pw3, 1);
pw8 = flip(pw2, 1);
N = 28; % sum(pw1(:));
PW = cat(3, pw1, pw2, pw3, pw4, pw5, pw6, pw7, pw8) / N;

% Filt each element
result = M.MapPage(@prewitt_filt, span, pw_id, PW, lookNum);

end


function result = prewitt_filt(cij, span, pw_id, PW, look_num)

result = zeros(size(cij));

span_2 = span.^2;

var_v = 1 / look_num;

for k = 1:8
    mask = pw_id == k;
    if ~any(mask, 'all')
        continue;
    end

    h = squeeze(PW(:,:,k));

    % Calculate stats by convolution
    z_mean_all = imfilter(span, h, 'replicate', 'same');
    z_sq_mean_all = imfilter(span_2, h, 'replicate', 'same');
    z_mean = z_mean_all(mask);
    z_sq_mean = z_sq_mean_all(mask);
    var_z = z_sq_mean - abs(z_mean).^2;

    % Calculate 'b'
    var_x = (var_z - z_mean.^2 * var_v) / (1 + var_v);
    b = (var_x + 1e-30) ./ (var_z + 1e-30);
    
    % Filt C_ij
    c_mean_all = imfilter(cij, h, 'replicate', 'same');
    c_mean = c_mean_all(mask);
    c_val = cij(mask);
    
    % set result
    result(mask) = c_mean + b .* (c_val - c_mean);
end

end

