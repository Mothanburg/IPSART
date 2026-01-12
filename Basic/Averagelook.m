function result = Averagelook(image, rowLook, colLook)

arguments
    image (:,:) double
    rowLook double
    colLook double
end

if rowLook == 1 && colLook == 1
    result = image;
    return;
end

h = fspecial("average", [rowLook colLook]);
result = imfilter(image, h);

end