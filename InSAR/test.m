[m,s]=StackCreating(dataMaster, [9582 1956], dataSlave, [9573 1956]);

mst = m.image;
slv = s.image;

[gcpm,gcps,gcpcorr]=CrossCorrRegistrating(mst, slv, FINE_WINDOW_SIZE=[128 128]);

out = Warping(slv, gcpm, gcps, gcpcorr, 2, "Sinc");

co = CoherenceComputing(mst,out,5);
figure; imagesc(co);

figure;
imagesc(db(abs(mst)));
hold on; scatter(gcpm(:,2), gcpm(:,1), '.', 'CData', gcpcorr);

figure;
imagesc(db(abs(slv)));
hold on; scatter(gcps(:,2), gcps(:,1), '.', 'CData', gcpcorr);