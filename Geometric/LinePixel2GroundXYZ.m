function xyz = LinePixel2GroundXYZ(line, pixel, centerXYZ, orbitCoeff, prf, samplingFreq, azStartTime, azEndTime, rangeDelay, options)

arguments
    line (1,1)
    pixel (1,1)
    centerXYZ (1,3)
    orbitCoeff (:,6)
    prf (1,1)
    samplingFreq (1,1)
    azStartTime (1,1)
    azEndTime (1,1)
    rangeDelay (1,1)
    options.MAX_ITER = 20
    options.TOLERANCE = 1e-16
end

az_time = Line2AzimuthTime(line, azStartTime, prf);
rg_time = Pixel2RangeTime(pixel, rangeDelay);

sat_state = OrbitInterp(orbitCoeff, azimuthTime, azStartTime, azEndTime);
sat_pos = sat_state(1:3);
sat_vel = sat_state(4:6);

xyz = centerXYZ;

for it = 1:options.MAX_ITER
    vec_inc = xyz - sat_pos;
    f_xyz = [
        -dot(sat_vel, vec_inc);
        -(dot(vec_inc, vec_inc) - (299792458 * rg_time)^2);
        -((xyz(1)^2 + xyz(2)^2)/ 6356752.31424518^2 + xyz(3)^2 / 6378137^2 - 1)
    ];
    jacobian_xyz = [
        sat_vel;
        2 * vec_inc;
        2 * xyz ./ [6356752.31424518 6356752.31424518 6378137].^2
    ];

    d_xyz = jacobian_xyz \ f_xyz;

    xyz = xyz + d_xyz;

    if all(abs(d_xyz)) < options.TOLERANCE
        break;
    end
end

end