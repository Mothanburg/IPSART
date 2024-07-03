function classes = SVM(features, sample, label)

model = fitcecoc(sample, label);

[height,width,~] = size(features);
classes = zeros(height, width);
len = height*width;
parfor i=1:len
    [m,n] = ind2sub([height width], i);
    f = squeeze(features(m,n,:))';
    classes(i) = model.predict(f);
end

end