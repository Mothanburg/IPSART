function result = AverageFilting(image, rowLook, colLook)

arguments
    image (:,:) double
    rowLook double
    colLook double
end

h = fspecial("average", [rowLook colLook]);
result = imfilter(image, h);

end