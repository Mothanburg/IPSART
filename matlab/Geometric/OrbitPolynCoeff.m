function varargout = OrbitPolynCoeff(orbitVector, order, startTime, endTime)

n_records = length(orbitVector);
if n_records < order - 1
    error("The order of polynomial is larger than amount of orbit records");
end

orbit_mat = zeros(n_records, 6);
x = zeros(n_records, 1);
for idx = 1:length(orbitVector)
    x(idx) = 4 * (orbitVector(idx).AzimuthTime - startTime) / (endTime - startTime) - 2;
    orbit_mat(idx, 1) = orbitVector(idx).PosX;
    orbit_mat(idx, 2) = orbitVector(idx).PosY;
    orbit_mat(idx, 3) = orbitVector(idx).PosZ;
    orbit_mat(idx, 4) = orbitVector(idx).VelX;
    orbit_mat(idx, 5) = orbitVector(idx).VelY;
    orbit_mat(idx, 6) = orbitVector(idx).VelZ;
end

A = ones(n_records, order + 1);
for n = 1:order
    A(:,n + 1) = x.^n;
end

coeff = pinv(A) * orbit_mat;

if nargout == 1
    varargout{1} = coeff;
elseif nargout == 2
    varargout{1} = coeff(:,1:3);
    varargout{2} = coeff(:,4:6);
else
    error("Invalid output arguments num");
end

end