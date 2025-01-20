function rangeTime = Pixel2RangeTime(pixel, rangeInitTime, samplingFreq)

rangeTime = rangeInitTime + (pixel - 1) / (2 * samplingFreq);

end