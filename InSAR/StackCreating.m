function [varargout] = StackCreating(master, masterPSPoint, slave, slavePSPoint)

arguments
    master
    masterPSPoint (1,2)
end

arguments (Repeating)
    slave
    slavePSPoint (1,2) {mustBeInteger}
end

[m_height,m_width] = size(master);

for i = 1:length(slave)
    [s_heights(i),s_widths(i)] = size(slave{i});

    offset = slavePSPoint{i} - masterPSPoint;
    pos_offset_hs(i) = offset(1);
    pos_offset_ws(i) = offset(2);

    neg_offset_hs(i) = s_heights(i) - pos_offset_hs(i) - m_height;
    neg_offset_ws(i) = s_widths(i) - pos_offset_ws(i) - m_width;
end

m_start_row = max([1, 1 - pos_offset_hs]);
m_end_row = min([m_height, m_height + neg_offset_hs]);
m_start_col = max([1, 1 - pos_offset_ws]);
m_end_col = min([m_width, m_width + neg_offset_ws]);

if m_start_row >= m_end_row || m_start_col >= m_end_col
    error("The images must have an overlapping area");
end

out.image = master(m_start_row:m_end_row,m_start_col:m_end_col);
out.row_skip = m_start_row - 1;
out.col_skip = m_start_col - 1;
varargout{1} = out;

for i = 1:length(slave)
    slave_img = slave{i};

    s_start_row = max([1, 1 + pos_offset_hs(i), 1 + pos_offset_hs(i) - min(pos_offset_hs)]);
    s_end_row = min([ ...
        s_heights(i), ...
        s_heights(i) - neg_offset_hs(i), ...
        s_heights(i) - neg_offset_hs(i) + min(neg_offset_hs) ...
        ]);

    s_start_col = max([1, 1 + pos_offset_ws(i), 1 + pos_offset_ws(i) - min(pos_offset_ws)]);
    s_end_col = min([ ...
        s_widths(i), ...
        s_widths(i) - neg_offset_ws(i), ...
        s_widths(i) - neg_offset_ws(i) + min(neg_offset_ws) ...
        ]);

    out.image = slave_img(s_start_row:s_end_row,s_start_col:s_end_col);
    out.row_skip = s_start_row - 1;
    out.col_skip = s_start_col - 1;
    varargout{i + 1} = out;
end


end