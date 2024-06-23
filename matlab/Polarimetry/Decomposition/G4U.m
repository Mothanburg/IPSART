function [Ps,Pd,Pv,Ph] = G4U(T3)

arguments
    T3 (:,:,3,3)
end

[height,width,~,~] = size(T3);
Ps = zeros(height, width);
Pd = zeros(height, width);
Pv = zeros(height, width);
Ph = zeros(height, width);

T3 = parallel.pool.Constant(shiftdim(T3, 2));
parfor j = 1:width
    for i = 1:height
        t0 = T3.Value(:,:,i,j);
        two_theta = atan(2 * real(t0(2,3)) / real(t0(2,2) - t0(3,3))) / 2;
        if isnan(two_theta)
            r = eye(3);
        else
            r = [1 0 0; 0 cos(two_theta) sin(two_theta); 0 -sin(two_theta) cos(two_theta)];
        end
        t = r * t0 * r';
        t11 = real(t(1,1));
        t22 = real(t(2,2));
        t33 = real(t(3,3));
        t12 = t(1,2);
        t13 = t(1,3);
        t23 = t(2,3);
        tp = t11 + t22 + t33;

        fh = 2 * abs(imag(t23));
        c1 = t11 - t22 + 7 * t33 / 8 + fh / 16;
        if c1 > 0
            co_ratio = 10 * log10((t11 + t22 - 2 * real(t12)) / (t11 + t22 + 2 * real(t12)));
            if co_ratio < -2
                fv = 15 * (2 * t33 - fh) / 8;
                if fv < 0
                    fh = 0;
                    fv = 15 * (2 * t33 - fh) / 8;
                end
                s = t11 - fv / 2;
                d = tp - fv - fh - s;
                c = t12 + t13 - fv / 6;
            elseif co_ratio > -2 && co_ratio < 2
                fv = 2 * (2 * t33 - fh);
                if fv < 0
                    fh = 0;
                    fv = 2 * (2 * t33 - fh);
                end
                s = t11 - fv / 2;
                d = tp - fv - fh - s;
                c = t12 + t13;
            else
                fv = 15 * (2 * t33 - fh) / 8;
                if fv < 0
                    fh = 0;
                    fv = 15 * (2 * t33 - fh) / 8;
                end
                s = t11 - fv / 2;
                d = tp - fv - fh - s;
                c = t12 + t13 + fv / 6;
            end

            if fv + fh > tp
                Ps(i,j) = 0;
                Pd(i,j) = 0;
                Ph(i,j) = fh;
                Pv(i,j) = tp - fh;
                continue;
            else
                c0 = 2 * t11 + fh - tp;
                if c0 > 0
                    fs = s + abs(c)^2 / s;
                    fd = d - abs(c)^2 / s;
                else
                    fs = s - abs(c)^2 / d;
                    fd = d + abs(c)^2 / d;
                end
            end
        else
            fv = 15 * (2 * t33 - fh) / 16;
            if fv < 0
                fh = 0;
                fv = 15 * (2 * t33 - fh) / 16;
            end
            s = t11;
            d = tp - fv - fh - s;
            c = t12 + t13;
            fs = s - abs(c)^2 / d;
            fd = d + abs(c)^2 / d;
        end

        if fs > 0 && fd > 0
            Ps(i,j) = fs;
            Pd(i,j) = fd;
            Ph(i,j) = fh;
            Pv(i,j) = fv;
        elseif fs > 0 && fd < 0
            Ps(i,j) = tp - fv - fh;
            Pd(i,j) = 0;
            Ph(i,j) = fh;
            Pv(i,j) = fv;
        elseif fs < 0 && fd > 0
            Ps(i,j) = 0;
            Pd(i,j) = tp - fh - fv;
            Ph(i,j) = fh;
            Pv(i,j) = fv;
        else
            error("未知错误");
        end

    end
end

end
