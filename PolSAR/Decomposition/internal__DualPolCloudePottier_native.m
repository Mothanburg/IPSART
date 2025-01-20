function [H,alpha,A] = internal__DualPolCloudePottier_native(M2)

if M2.Dtype == "double"
    [H,alpha,A] = clib.gars.CloudePottier2d(M2.m11, M2.m22, M2.m12_r, M2.m12_i);
else
    [H,alpha,A] = clib.gars.CloudePottier2f(M2.m11, M2.m22, M2.m12_r, M2.m12_i);
end

end