function C3 = fn_T3_to_C3(T3)

arguments
    T3 (:,:,3,3)
end

T3 = shiftdim(T3, 2);
M = [1 0 1; 1 0 -1; 0 sqrt(2) 0] / sqrt(2);

C3 = pagemtimes(pagemtimes(M', T3), M);
C3 = shiftdim(C3, 2);

end