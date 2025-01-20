function azimuthTime = Line2AzimuthTime(line, azimuthInitTime, prf)

azimuthTime = (line - 1) / prf + azimuthInitTime;

end
