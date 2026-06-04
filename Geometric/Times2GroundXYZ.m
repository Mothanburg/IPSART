function xyz = Times2GroundXYZ( ...
    azimuthTime, ...
    rangeTime, ...
    orbitPlyn, ...
    centerXYZ, ...
    ecefMajorAxis, ...
    ecefMinorAxis, ...
    options)

arguments
    azimuthTime (1,1) {mustBeInteger,mustBeGreaterThan(azimuthTime, 0)}
    rangeTime (1,1) {mustBeInteger,mustBeGreaterThan(rangeTime, 0)}
    orbitPlyn OrbitPolynomial
    centerXYZ (1,3) double
    ecefMajorAxis (1,1) double
    ecefMinorAxis (1,1) double
    options.MAX_ITER = 20
    options.TOLERANCE = 1e-16
end

sat_state = orbitPlyn.GetSataliteState(azimuthTime);
sat_pos = sat_state(1:3);
sat_vel = sat_state(4:6);

xyz = centerXYZ;

for it = 1:options.MAX_ITER
    vec_inc = xyz - sat_pos;
    f_xyz = [
        -dot(sat_vel, vec_inc);
        -(dot(vec_inc, vec_inc) - (299792458 * rangeTime)^2);
        -((xyz(1)^2 + xyz(2)^2)/ ecefMinorAxis^2 + xyz(3)^2 / ecefMajorAxis - 1)
        ];
    jacobian_xyz = [
        sat_vel;
        2 * vec_inc;
        2 * xyz ./ [ecefMinorAxis ecefMinorAxis ecefMajorAxis].^2
        ];

    d_xyz = jacobian_xyz \ f_xyz;

    xyz = xyz + d_xyz;

    if all(abs(d_xyz)) < options.TOLERANCE
        break;
    end
end

end 