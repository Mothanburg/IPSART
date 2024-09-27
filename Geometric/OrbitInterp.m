function result = OrbitInterp(orbitCoeff, azimuthTime, startTime, endTime)

[order_p1,~] = size(orbitCoeff);
azimuthTime = 4 * (azimuthTime - startTime) / (endTime - startTime) - 2;
A = ones(length(azimuthTime), order_p1);
for n = 2:order_p1
    A(:,n) = azimuthTime.^(n - 1);
end

result = A * orbitCoeff;

end