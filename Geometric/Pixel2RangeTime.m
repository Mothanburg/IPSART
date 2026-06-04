% Pixel2RangeTime - 将图像的列数转换为对应的距离向采集时间
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function rangeTime = Pixel2RangeTime(pixel, rangeInitTime, samplingFreq)

rangeTime = rangeInitTime + (pixel - 1) / (2 * samplingFreq);

end