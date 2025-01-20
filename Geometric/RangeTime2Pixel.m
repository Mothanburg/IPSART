function pixel = RangeTime2Pixel(rangeTime, rangeInitTime, samplingFreq)

pixel = 2 * (rangeTime - rangeInitTime) * samplingFreq + 1;

end