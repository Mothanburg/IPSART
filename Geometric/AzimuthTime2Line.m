function line = AzimuthTime2Line(azimuthTime, azimuthInitTime, prf)

line = (azimuthTime - azimuthInitTime) * prf + 1;

end
