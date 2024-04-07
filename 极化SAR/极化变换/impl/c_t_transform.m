function result = c_t_transform(input, t2c, dataSpec)

[height,width,dim,~] = size(input);
switch dataSpec
    case {'hhhv', 'vvvh'}
        u = [sqrt(2), 0; 0, 1+1i]/sqrt(2);
    case 'hhvv'
        u = [1, 1; 1, -1];
    case 'quad'
        u = [1,0,1;1,0,-1;0,sqrt(2),0]/sqrt(2);
    otherwise
        error('错误的数据类型');
end

if t2c
    u = u';
end

len = height*width;
result = zeros(height, width, dim, dim);
for i=1:len
    [m,n] = ind2sub([height width], i);
    m_in = squeeze(input(m,n,:,:));
    m_out = u*m_in*u';
    result(m,n,:,:) = m_out;
end

end