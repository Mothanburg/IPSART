function result = Averagelook(image, rowLook, colLook)

arguments
    image (:,:) double
    rowLook double
    colLook double
end

h = fspecial("average", [rowLook colLook]);
result = imfilter(image, h);

end