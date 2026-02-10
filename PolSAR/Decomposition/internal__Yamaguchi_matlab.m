function [Ps,Pd,Pv,Ph] = internal__Yamaguchi_matlab(T3)

height = T3.Height;
width = T3.Width;

Ps = zeros(height, width, T3.Dtype);
Pd = zeros(height, width, T3.Dtype);
Pv = zeros(height, width, T3.Dtype);
Ph = zeros(height, width, T3.Dtype);

parfor j = 1:width
    for i = 1:height
        t0 = T3.MatAt(i, j);

        % matrix rotation
        theta = atan2(2 * real(t0(2,3)), t0(2,2) - t0(3,3)) / 4;
        if theta < -pi / 4
            theta = theta + pi / 2;
        elseif theta > pi / 4
            theta = theta - pi / 2;
        elseif isnan(theta) % theta is nan
            theta = 0;
        end
        
        R = [1 0 0; 0 cos(2 * theta) sin(2 * theta); 0 -sin(2 * theta) cos(2 * theta)];
        t = R * t0 * R';

        tp = sum(real(diag(t)));

        % helix component
        fh = 2 * abs(imag(t(2,3)));

        % dominant mechanism
        if t(1,1) > t(2,2) - fh / 2 % surface
            % volume scattering form
            ratio = 10 * log10((t(1,1) + t(1,2) - 2 * real(t(1,2))) / (t(1,1) + t(1,2) + 2 * real(t(1,2))));
            if ratio < -2
                fv = 15 * (t(3,3) / 4 - fh / 8);
                if fv < 0 % remove helix if volume is too big
                    fh = 0;
                    fv = 15 * (t(3,3) / 4 - fh / 8);
                end
                S = t(1,1) - fv / 2;
                D = t(2,2) - 7 * fv / 30 - fh / 2;
                C = t(1,2) - fv / 6;
            elseif ratio < 2
                fv = 4 * t(3,3) - 2 * fh;
                if fv < 0 % remove helix if volume is too big
                    fh = 0;
                    fv = 4 * t(3,3) - 2 * fh;
                end
                S = t(1,1) - fv / 2;
                D = t(2,2) - t(3,3);
                C = t(1,2);
            else
                fv = 15 * (t(3,3) / 4 - fh / 8);
                if fv < 0 % remove helix if volume is too big
                    fh = 0;
                    fv = 15 * (t(3,3) / 4 - fh / 8);
                end
                S = real(t(1,1)) - fv / 2;
                D = real(t(2,2)) - 7 * fv / 30 - fh / 2;
                C = t(1,2) + fv / 6;
            end

            % caculate ps and pd
            if fv + fh > tp % no ps and pd
                fs = 0;
                fd = 0;
                fv = tp - fh;

                % to end
                Ps(i,j) = fs;
                Pd(i,j) = fd;
                Ph(i,j) = fh;
                Pv(i,j) = fv;
                continue;
            else
                C0 = t(1,1) - t(2,2) - t(3,3) + fh;
                if C0 > 0
                    fs = S + abs(C)^2 / S;
                    fd = D - abs(C)^2 / S;
                else % goto double bounce
                    fs = S - abs(C)^2 / D;
                    fd = D + abs(C)^2 / D;
                end
            end

        else % double bounce
            fv = 15 * (t(3,3) - fh / 2) / 8;
            if fv < 0 % remove helix if volume is too big
                fh = 0;
                fv = 15 * (t(3,3) - fh / 2) / 8;
            end
            S = real(t(1,1));
            D = real(t(2,2)) - 7 * fv / 15 - fh / 2;
            C = t(1,2);

            fs = S - abs(C)^2 / D;
            fd = D + abs(C)^2 / D;
        end

        if fs > 0 && fd < 0
            fd = 0;
            fs = tp - fv - fh;
        elseif fs < 0 && fd > 0
            fs = 0;
            fd = tp - fv - fh;
        elseif fs < 0 && fd < 0 % impossible
            error("Unknown error");
        end

        Ps(i,j) = real(fs);
        Pd(i,j) = real(fd);
        Ph(i,j) = real(fh);
        Pv(i,j) = real(fv);
    end
end

end
