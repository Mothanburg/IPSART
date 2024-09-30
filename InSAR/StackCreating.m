function [mstOut,varargout] = StackCreating(master, slaves, slaveOffsets)

arguments
    master
end

arguments (Repeating)
    slaves
    slaveOffsets (1,2) {mustBeInteger}
end

[mlines,mpixels] = size(master);

for i = 1:length(slaves)
    [sliness(i),spixelss(i)] = size(slaves{i});

    offset = slaveOffsets{i};
    offsetls(i) = offset(1);
    offsetps(i) = offset(2);

    rev_offsetls(i) = sliness(i) - offsetls(i) - mlines;
    rev_offsetps(i) = spixelss(i) - offsetps(i) - mpixels;
end

mststartl = max([1, 1 - offsetls]);
mstendl = min([mlines, mlines + rev_offsetls]);
mststartp = max([1, 1 - offsetps]);
mstendp = min([mpixels, mpixels + rev_offsetps]);

if mststartl >= mstendl || mststartp >= mstendp
    error("The images must have an overlapping area");
end

mstOut = master(mststartl:mstendl,mststartp:mstendp);

for i = 1:length(slaves)
    slave = slaves{i};

    slvstartl = max([1, 1 + offsetls(i), 1 + offsetls(i) - min(offsetls)]);
    slvendl = min([ ...
        sliness(i), ...
        sliness(i) - rev_offsetls(i), ...
        sliness(i) - rev_offsetls(i) + min(rev_offsetls) ...
        ]);

    slvstartp = max([1, 1 + offsetps(i), 1 + offsetps(i) - min(offsetps)]);
    slvendp = min([ ...
        spixelss(i), ...
        spixelss(i) - rev_offsetps(i), ...
        spixelss(i) - rev_offsetps(i) + min(rev_offsetps) ...
        ]);

    varargout{i} = slave(slvstartl:slvendl,slvstartp:slvendp);
end


end