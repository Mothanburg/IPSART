function rangeTime = Pixel2RangeTime(pixel, rangeDelay, samplingFreq)

    rangeTime = rangeDelay + (pixel - 1) / (2 * samplingFreq);

end