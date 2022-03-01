function CoM = AdaptiveThresh(img, gsig)

filt = (fspecial('gaussian', [gsig gsig]));
out = conv2(img, filt, 'same');

im = (out - min(out,[],'all')) / (max(out,[],'all') - min(out,[],'all')) * 255;
im = uint8(im);

thr = adaptthresh(im, 'NeighborhoodSize', gsig, 'Statistic', 'gaussian') ;
bin = imbinarize(im, thr);

th2 = bwareaopen(bin, 10);

th3 = bwmorph(th2, 'open');
th3 = bwareaopen(th3, 10);

labeled = bwlabel(th3);
blobCOM = regionprops(labeled, 'Centroid');
CoM = cell2mat(arrayfun(@(x) x.Centroid, blobCOM, 'UniformOutput', false));
CoM = round(CoM);