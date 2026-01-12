function [azimuthTime,rangeTime] = ...
    GroundXYZ2Times(groundPos, orbitPlyn, options)

arguments
    groundPos (1,3)
    orbitPlyn OrbitPolynomial
    options.MAX_ITER = 20
    options.TOLERANCE = 1e-16
end

azimuthTime = mean(orbitPlyn.TimeRange);
for it = 1:options.MAX_ITER
    sat_state = orbitPlyn.GetSataliteState(azimuthTime);
    sat_pos = sat_state(1:3);
    sat_v = sat_state(4:6);

    % acceleration of the satallite, need a derivation: a = dv / dt 
    % the degree of polynomial will reduce 1
    T = orbitPlyn.TimeMatrix(azimuthTime, orbitPlyn.Degree - 1);
    sat_acc = 4 * T * orbitPlyn.Coeff(2:end,4:6) / diff(orbitPlyn.TimeRange);

    vec_inc = groundPos - sat_pos;
    dt = -dot(sat_v, vec_inc) / (dot(sat_acc, vec_inc) - dot(sat_v, sat_v));
    azimuthTime = azimuthTime + dt;

    if abs(dt) <= options.TOLERANCE
        break;
    end
end

rangeTime = norm(vec_inc) / physconst("LightSpeed");

end
