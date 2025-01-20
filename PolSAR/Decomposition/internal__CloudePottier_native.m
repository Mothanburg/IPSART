function [H,alpha,A] = internal__CloudePottier_native(M3)

if M3.Dtype == "double"
    [H,alpha,A] = clib.gars.CloudePottier3d(M3.m11, M3.m22, M3.m33, ...
        M3.m12_r, M3.m13_r, M3.m23_r, M3.m12_i, M3.m13_i, M3.m23_i);
else
    [H,alpha,A] = clib.gars.CloudePottier3f(M3.m11, M3.m22, M3.m33, ...
        M3.m12_r, M3.m13_r, M3.m23_r, M3.m12_i, M3.m13_i, M3.m23_i);
end

end