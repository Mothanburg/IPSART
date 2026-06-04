% Line2AzimuthTime - 将图像的行数转换为对应的方位向采集时间
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function azimuthTime = Line2AzimuthTime(line, azimuthInitTime, prf)

azimuthTime = (line - 1) / prf + azimuthInitTime;

end
