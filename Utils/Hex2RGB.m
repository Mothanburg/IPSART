% Hex2RGB - 将HEX颜色字符串转换为MATLAB RGB数值
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function rgb = Hex2RGB(hex)
    hex = char(upper(hex));
    exchange_list='0123456789ABCDEF#';
    rgb = zeros(1,3);
    for i = 1:3
        tempCoe1=find(exchange_list==hex(i*2))-1;
        tempCoe2=find(exchange_list==hex(i*2+1))-1;
        rgb(i)=16*tempCoe1+tempCoe2;
    end
    rgb = rgb / 255;
end