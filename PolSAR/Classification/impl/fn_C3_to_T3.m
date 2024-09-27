function T3 = fn_C3_to_T3(C3)

arguments
    C3 (:,:,3,3)
end

C3 = shiftdim(C3, 2);
M = [1 0 1; 1 0 -1; 0 sqrt(2) 0] / sqrt(2);

T3 = pagemtimes(pagemtimes(M, C3), M');
T3 = shiftdim(T3, 2);

end