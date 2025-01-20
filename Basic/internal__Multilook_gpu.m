function result = internal__Multilook_gpu( ...
    image, in_rows, in_cols, row_look, col_look)

out_rows = fix(in_rows / row_look);
out_cols = fix(in_cols / col_look);

if ~isreal(image)
    if class(image) == "double"
        [errno1,rpart] = clib.gars.Multilookd(real(image), row_look, ...
            col_look, out_rows, out_cols);
        [errno2,ipart] = clib.gars.Multilookd(imag(image), row_look, ...
            col_look, out_rows, out_cols);
    else
        [errno1,rpart] = clib.gars.Multilookf(real(image), row_look, ...
            col_look, out_rows, out_cols);
        [errno2,ipart] = clib.gars.Multilookf(imag(image), row_look, ...
            col_look, out_rows, out_cols);
    end

    if errno1 ~= 0 || errno2 ~= 0
        error("error number %d and %d is returned", errno1, errno2);
    end

    result = rpart + 1i * ipart;
else
    if class(image) == "double"
        [errno,result] = clib.gars.Multilookd(image, row_look, col_look, ...
            out_rows, out_cols);
    else
        [errno,result] = clib.gars.Multilookf(image, row_look, col_look, ...
            out_rows, out_cols);
    end

    if errno ~= 0
        error("error number %d is returned", errno);
    end
end

end