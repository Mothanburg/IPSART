% AzimuthTime2Line - 将方位向采集转换为图像上对应的行数
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function line = AzimuthTime2Line(azimuthTime, azimuthInitTime, prf)

line = (azimuthTime - azimuthInitTime) * prf + 1;

end
