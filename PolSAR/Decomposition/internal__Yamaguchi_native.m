function [Ps,Pd,Pv,Ph] = internal__Yamaguchi_native(C3)

if C3.Dtype == "double"
    [Ps,Pd,Pv,Ph] = clib.gars.Yamaguchid(C3.m11, C3.m22, C3.m33, C3.m12_r, ...
        C3.m13_r, C3.m23_r, C3.m12_i, C3.m13_i, C3.m23_i);
else
    [Ps,Pd,Pv,Ph] = clib.gars.Yamaguchif(C3.m11, C3.m22, C3.m33, C3.m12_r, ...
        C3.m13_r, C3.m23_r, C3.m12_i, C3.m13_i, C3.m23_i);
end

end