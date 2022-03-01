function [SourceCell]=extractROIs3(img,Opts)
sources=[];

if ~isfield(Opts,'radius')
radius = 5;
else
radius = Opts.radius;
end

if ~isfield(Opts,'SizeThreshold')
SizeThreshold = 250;
else
SizeThreshold = Opts.SizeThreshold;
end

if ~isfield(Opts,'maxSourcesPerPlane')
MaxSources = 250;
else
MaxSources = Opts.maxSourcesPerPlane;
end

if ~isfield(Opts,'height')
height=1:512;
else
height = Opts.height;
end

if ~isfield(Opts,'width')
width=100:400;
else
width=Opts.width;
end

if ~isfield(Opts,'distThreshold')
distThreshold=0;
preventOverlap=0;
else
preventOverlap=1;
distThreshold=Opts.distThreshold;;

end

nd=size(img,3);



thresh=[];

% preprocess
% filt=(fspecial('gaussian',[30 30],2));

for j = 1:nd  % for each depth
    z=zeros(512,512);
    temp=squeeze(img(:,:,j));
%     temp2=conv2(temp,filt,'same');
    temp2 = temp;
    z(height,width)=temp2(height,width);
    Data(j).raw=z;
    Data(j).thresh=mean(temp(:))+std(temp(:));
end;

imageSizeX = 512;
imageSizeY = 512;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
for k = 1:nd;
    clear OC;
    tic
    disp(['Segmenting Depth ' num2str(k)]);
    
    [output]  =AdaptiveThresh(Data(k).raw, 5);
    if(isempty(output))
        output = ones(4,1);
    end
    OC(:,1)=output(1:2:end);
    OC(:,2)=output(2:2:end);
    
    for j = 1:length(OC)
        circlePixels = (rowsInImage  - OC(j,2)).^2 + (columnsInImage - OC(j,1)).^2 <= radius.^2;
        Flo=circlePixels.*Data(k).raw;
        Flo=Flo(:);
        Flo(Flo==0)=[];
        OC(j,3)=mean(Flo);
    end
    
    
    OC(OC(:,3)<Data(k).thresh,:)=[];
    
    if Opts.brightestFirst
    OC=sortrows(OC,3);
    OC=flipud(OC);    %brightest first
    else
        OC=OC(randperm(size(OC,1)),:);
    end
    
    if preventOverlap
        distMatrix =[];
        for a=1:size(OC,1)
            for b = 1:size(OC,1);
                distMatrix(a,b) = sqrt(mean((OC(a,1:2) - OC(b,1:2)).^2));
                if a ==b
                    distMatrix(a,b)=nan;
                end
            end
        end
        
        keepList=[];
        for j = 1:size(OC,1)
            dists = distMatrix(j,:);
            if  ~any(dists(keepList)<distThreshold)
                keepList = [keepList j];
            end
        end
        OC=OC(keepList,:);
    end
    
    
    if size(OC,1)>MaxSources;
        try
            Data(k).OC=OC(1:MaxSources,:);
        catch
            Data(k).OC=OC;
        end
    else
        Data(k).OC=OC;
    end
    
   sources = zeros(512,512,size(Data(k).OC,1));
   SE=strel('disk',radius,4);
      
   for n = 1:size(sources,3);
         sources(Data(k).OC(n,2),Data(k).OC(n,1),n)=1;
         sources(:,:,n)=imdilate(sources(:,:,n),SE);
   end
    
   SourceCell{k}=sources;     
    
    toc
end