function azimuthTime = Line2AzimuthTime(line, startTime, prf)

  azimuthTime = (line - 1) / prf + startTime;

end
