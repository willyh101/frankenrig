%% Click to get targets
clear imgClick
imgClick{1}=ImageData;
imgClick{2}=ImageData;
imgClick{3}=ImageData;
%%
nImg = numel(imgClick);
% zplanes = hSI.hStackManager.arbitraryZs;
zplanes = [0];
nplanes = numel(zplanes);
button=1;
locs = [];

f = figure(1111);
for i=1:nImg
    clf
    image(imgClick{i})
    z = zplanes(i);
    title(['Plane ' num2str(z) '. Press enter to continue.'])
    [x,y,button] = ginput;
    if numel(button) == 0
        locs{i} = {};
    else
        locs{i}(:,1:2) = round([x y]);
        locs{i}(:,3) = ones(numel(x),1)*z;
    end 
end
close(f)


%%make sources

radius = 5;
SE=strel('disk',radius,4);

for k=1:nplanes
    sources_temp = zeros(512,512,size(locs{k},1));
    for n = 1:size(sources_temp,3)
       sources_temp(locs{k}(n,2),locs{k}(n,1),n)=1;
       sources_temp(:,:,n)=imdilate(sources_temp(:,:,n),SE);
    end
    all_sources{k} = sources_temp;
end


%%send targets to ScanImage 

hSI.hIntegrationRoiManager.roiGroup.clear()
imagingScanfield = hSI.hRoiManager.currentRoiGroup.rois(1).scanfields(1);

correctZ =1; %Hacky correct Z errors added 9/30/19 by Ian
HoloOffsets;

number = 0;
for k = 1:nplanes
    sources = all_sources{k};
    for i = 1:size(sources,3)
       mask = sources(:,:,i);
       intsf = scanimage.mroi.scanfield.fields.IntegrationField.createFromMask(imagingScanfield,mask);
       intsf.threshold = 100;
       introi = scanimage.mroi.Roi();
       introi.discretePlaneMode=1;
       introi.add(zplanes(k), intsf);
       introi.name = ['ROI ' num2str(number+1) ];%' Depth ' num2str(zDepth(n))];
       hSI.hIntegrationRoiManager.roiGroup.add(introi);
       number=number+1;
    end
end

disp(['Added ' num2str(number) ' sources to integration']);

holoRequests = saveAllToHoloRequest(yoffset, xoffset,zMap); % y,x offsets (-4, -2.8)
selectAllROIs;

disp('sent ROIs to the cloud');