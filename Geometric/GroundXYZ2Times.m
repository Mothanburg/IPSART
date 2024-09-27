function varargout = ...
    GroundXYZ2Times(groundPos, orbitCoeff, startTime, endTime, options)

arguments
    groundPos (1,3)
    orbitCoeff (:,6)
    startTime (1,1)
    endTime (1,1)
    options.MAX_ITER = 20
    options.TOLERANCE = 1e-16
end

order = numel(orbitCoeff) / 6 - 1;

azimuthTime = (startTime + endTime) / 2;
for it = 1:options.MAX_ITER
    sat_state = OrbitInterp(orbitCoeff, azimuthTime, startTime, endTime);
    sat_pos = sat_state(1:3);
    sat_vel = sat_state(4:6);
    sat_acc = 4 * OrbitInterp( ...
        orbitCoeff(2:end,4:6) .* (1:order)', azimuthTime, startTime, endTime) / ...
        (endTime - startTime);

    vec_inc = groundPos - sat_pos;
    dt = -dot(sat_vel, vec_inc) / (dot(sat_acc, vec_inc) - dot(sat_vel, sat_vel));
    azimuthTime = azimuthTime + dt;

    if abs(dt) <= options.TOLERANCE
        break;
    end
end

if nargout == 1
    varargout{1} = azimuthTime;
elseif nargout == 2
    varargout{1} = azimuthTime;
    varargout{2} = norm(vec_inc) / 299792458;
else
    error("Invalid count of output arguments");
end

end