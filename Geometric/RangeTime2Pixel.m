function pixel = RangeTime2Pixel(rangeTime, rangeDelay, samplingFreq)

    pixel = 2 * (rangeTime - rangeDelay) * samplingFreq + 1;

end