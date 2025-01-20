function result = internal__RefinedLeeFilting_gpu(M, lookNum)

if isa(M, "PolM2")
    result = rlf_2x2(M, lookNum);
else
    result = rlf_3x3(M, lookNum);
end

end


function result = rlf_2x2(M, lookNum)

if M.Dtype == "double"
    [errno, m11, m22, m12r, m12i] = ...
        clib.gars.RefinedLeeFilting2d(lookNum, M.m11, M.m22, M.m12_r, M.m12_i);
else
    [errno, m11, m22, m12r, m12i] = ...
        clib.gars.RefinedLeeFilting2f(lookNum, M.m11, M.m22, M.m12_r, M.m12_i);
end

if errno ~= 0
    error("error number %d is returned", errno);
end

if isa(M, "PolC2")
    result = PolC2(m11, m22, m12r, m12i, M.PolType);
elseif isa(M, "PolT2")
    result = PolT2(m11, m22, m12r, m12i);
end

end


function result = rlf_3x3(M, lookNum)

if M.Dtype == "double"
    [errno, m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = ...
        clib.gars.RefinedLeeFilting3d(lookNum, M.m11, M.m22, M.m33, ...
        M.m12_r, M.m13_r, M.m23_r, M.m12_i, M.m13_i, M.m23_i);
else
    [errno, m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i] = ...
        clib.gars.RefinedLeeFilting3f(lookNum, M.m11, M.m22, M.m33, ...
        M.m12_r, M.m13_r, M.m23_r, M.m12_i, M.m13_i, M.m23_i);
end

if errno ~= 0
    error("error number %d is returned", errno);
end

if isa(M, "PolC3")
    result = PolC3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
elseif isa(M, "PolT3")
    result = PolT3(m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i);
end

end