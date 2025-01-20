function [Ps,Pd,Pv,Ph] = internal__G4U_native(T3)

if T3.Dtype == "double"
    [errno,Ps,Pd,Pv,Ph] = clib.gars.G4Ud(T3.m11, T3.m22, T3.m33, T3.m12_r, ...
        T3.m13_r, T3.m23_r, T3.m12_i, T3.m13_i, T3.m23_i);
else
    [errno,Ps,Pd,Pv,Ph] = clib.gars.G4Uf(T3.m11, T3.m22, T3.m33, T3.m12_r, ...
        T3.m13_r, T3.m23_r, T3.m12_i, T3.m13_i, T3.m23_i);
end

if errno ~= 0
    warning("There may be abnormal values, please check");
end

end