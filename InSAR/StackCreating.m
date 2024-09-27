function [outMst,outSlv] = StackCreating(master, slave, slaveOffset, options)

arguments
    master
    slave
    slaveOffset
    options.OUTPUT_AREA = "overlap"
end

[mlines,mpixels] = size(master);
[slines,spixels] = size(slave);

switch lower(options.OUTPUT_AREA)
    case "overlap"
        outMst = master(...
            max(1, 1 - slaveOffset(1)):min(mlines, mlines - slaveOffset(1)), ...
            max(1, 1 - slaveOffset(2)):min(mpixels, mpixels - slaveOffset(2)) ...
            );
        outSlv = slave(...
            max(1, 1 + slaveOffset(1)):min(slines, slines + slaveOffset(1)), ...
            max(1, 1 + slaveOffset(2)):min(spixels, spixels + slaveOffset(2)) ...
            );
    case "master"
        outMst = master;
        outSlv = nan(mlines, mpixels);

        out_l0 = max(1, 1 - slaveOffset(1));
        out_p0 = max(1, 1 - slaveOffset(2));

        slv_l0 = max(1, 1 + slaveOffset(1));
        slv_p0 = max(1, 1 + slaveOffset(2));

        if slaveOffset(1) < 0
            cp_len_l = min(slines, mlines + slaveOffset(1));
        else
            cp_len_l = min(mlines, slines - slaveOffset(1));
        end
        
        if slaveOffset(2) < 0
            cp_len_p = min(spixels, mpixels + slaveOffset(2));
        else
            cp_len_p = min(mpixels, spixels - slaveOffset(2));
        end

        out_lines = out_l0 - 1 + (1:cp_len_l);
        out_pixels = out_p0 - 1 + (1:cp_len_p);
        slv_lines = slv_l0 - 1 + (1:cp_len_l);
        slv_pixels = slv_p0 - 1 + (1:cp_len_p);
        outSlv(out_lines,out_pixels) = slave(slv_lines,slv_pixels);
    otherwise
        error("Invalid options: %s", options.OUTPUT_AREA);
end

end