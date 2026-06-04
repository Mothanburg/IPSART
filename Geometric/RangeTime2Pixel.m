% RangeTime2Pixel - 将距离向采集时间转换为对应的图像列数
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function pixel = RangeTime2Pixel(rangeTime, rangeInitTime, samplingFreq)

pixel = 2 * (rangeTime - rangeInitTime) * samplingFreq + 1;

end