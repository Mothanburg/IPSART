function [Ps,Pd,Pv,Ph] = Yamaguchi(C3)

arguments
    C3 PolC3
end

global GARS_CONFIG

if GARS_CONFIG.CAPABILITY > 0
    try
        [Ps,Pd,Pv,Ph] = clib.gars.Yamaguchi(C3.m11, C3.m22, C3.m33, C3.m12_r, C3.m13_r, ...
            C3.m23_r, C3.m12_i, C3.m13_i, C3.m23_i);
        return
    catch e
        warning(e.identifier, "An error occurred when calling library, fallback to matlab.\n" + ...
            "    Error message: %s", e.message);
    end
end

[Ps,Pd,Pv,Ph] = Ymg_matlab(C3);

end

function [Ps,Pd,Pv,Ph] = Ymg_matlab(C3)

height = C3.Height;
width = C3.Width;

Ps = zeros(height, width, C3.Dtype);
Pd = zeros(height, width, C3.Dtype);
Pv = zeros(height, width, C3.Dtype);
Ph = zeros(height, width, C3.Dtype);

C3 = parallel.pool.Constant(C3);
parfor j = 1:width
    for i = 1:height
        c = C3.Value.getMatAt(i, j);

        % 计算螺旋体方向
        if imag(c(1,2) + c(2,3)) > 0
            Ch = [1 1i*sqrt(2) -1; -1i*sqrt(2) 2 1i*sqrt(2); -1 -1i*sqrt(2) 1] / 4; % 右旋
        else
            Ch = [1 -1i*sqrt(2) -1; 1i*sqrt(2) 2 -1i*sqrt(2); -1 1i*sqrt(2) 1] / 4; % 左旋
        end
        fh = 2 * abs(imag(c(1,2) + c(2,3)));

        % 计算体散射形式
        co_ratio = 10 * log10(real(c(3,3)) / real(c(1,1)));
        if co_ratio < -2
            Cv = [8 0 2; 0 4 0; 2 0 3] / 15;
            fv = 15 * (real(c(2,2)) - fh / 2) / 4;
        elseif co_ratio < 2
            Cv = [3 0 1; 0 2 0; 1 0 3] / 8;
            fv = 4 * (real(c(2,2)) - fh / 2);
        else
            Cv = [3 0 2; 0 4 0; 2 0 8] / 15;
            fv = 15 * (real(c(2,2)) - fh / 2) / 4;
        end

        % 仅在交叉通道足够小时减去体散射和螺旋散射
        if real(c(2,2)) < real(c(1,1)) && real(c(2,2)) < real(c(3,3))
            c = c - fh * Ch - fv * Cv;
        else
            fh = 0;
            fv = 0;
        end

        if real(c(1,3)) > 0
            a = -1;
            fd = real((c(3,3) * c(1,1) - c(1,3) * c(3,1)) / (c(3,3) + c(1,1) + c(1,3) + c(3,1)));
            fs = real(c(3,3)) - fd;
            b = (c(1,3) + fd) / fs;
        else
            b = 1;
            fs = real((c(1,3) * c(3,1) - c(3,3) * c(1,1)) / (c(1,3) + c(3,1) - c(1,1) - c(3,3)));
            fd = real(c(3,3)) - fs;
            a = (c(1,3) - fs) / fd;
        end

        Ps(i,j) = abs(fs) * (1 + b * conj(b));
        Pd(i,j) = abs(fd) * (1 + a * conj(a));
        Ph(i,j) = abs(fh);
        Pv(i,j) = abs(fv);
    end
end

end
