clear; close all; clc;
%%
offsetNum = 1;
%%
tiff_info = imfinfo(strcat('_',num2str(offsetNum),'\_',num2str(offsetNum),'_MMStack_Default.ome.tif')); % return tiff structure, one element per image
tiff_stack = imread(strcat('_',num2str(offsetNum),'\_',num2str(offsetNum),'_MMStack_Default.ome.tif'), 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(strcat('_',num2str(offsetNum),'\_',num2str(offsetNum),'_MMStack_Default.ome.tif'), ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

offsetTif = double(tiff_stack);
offset = mean(offsetTif,'all');
%%
fileNum = 11;

fovSize = 121;
centerLoc = [507, 505];

Radius = 45;

Dig2Ph = 0.29;
%%
tiff_info = imfinfo(strcat('_',num2str(fileNum),'\_',num2str(fileNum),'_MMStack_Default.ome.tif')); % return tiff structure, one element per image
tiff_stack = imread(strcat('_',num2str(fileNum),'\_',num2str(fileNum),'_MMStack_Default.ome.tif'), 10) ; % read in first image
%concatenate each successive tiff to tiff_stack
% for ii = 2 : min(size(tiff_info, 1),100)
%     temp_tiff = imread(strcat('_',num2str(fileNum),'\_',num2str(fileNum),'_MMStack_Default.ome.tif'), ii);
%     tiff_stack = cat(3 , tiff_stack, temp_tiff);
% end
%%
dataTif = double(tiff_stack);
avgDataTif = mean(dataTif,3)-offset;
avgDataTif = avgDataTif*Dig2Ph;
%% Cropping and computing
cropData = avgDataTif(centerLoc(1)-(fovSize-1)/2: centerLoc(1)+(fovSize-1)/2,...
    centerLoc(2)-(fovSize-1)/2:centerLoc(2)+(fovSize-1)/2);
centerLoc_new = (fovSize-1)/2 + 1;
[x,y]=meshgrid(-(fovSize-1)/2:(fovSize-1)/2);
radiusData = sqrt(x.^2+y.^2);
radiusNorm = radiusData./46;
%%
% r_tmp = reshape(radiusNorm,1,[]);
% i_tmp = reshape(cropData,1,[]);
% [R_tmp, ind]=sort(r_tmp);
% I_tmp = i_tmp(ind);
% plot(R_tmp,I_tmp);
%%
Radius_sample = 0:0.05:1.5;
radius_ind =  (Radius_sample(2:end)+Radius_sample(1:end-1))/2;
I_mean = zeros(size(radius_ind));
I_std =  zeros(size(radius_ind));
I_med =  zeros(size(radius_ind));
I_90 = zeros(size(radius_ind));
I_10 =  zeros(size(radius_ind));
for i = 1:length(Radius_sample)-1
    DataPc = cropData(radiusNorm>=Radius_sample(i) & radiusNorm<Radius_sample(i+1));
    I_mean(i) = mean(DataPc,'all');
    I_std(i) =  std(DataPc,0,'all');
    I_med(i) = prctile(DataPc,50,'all');
    I_90(i) = prctile(DataPc,90,'all');
    I_10(i) =  prctile(DataPc,10,'all');
end
%%
figure('Position',[100,300,300,450]); tiledlayout('flow'); hold on;
nexttile; imagesc(cropData); colormap gray; hcb = colorbar; axis image; box off; axis off; 
line([100 100+17.0940],[110 110],'Color','w','LineWidth',4);
text(102,115,'1um','Color','w','FontWeight','bold');

set(get(hcb,'Title'),'String','photons','fontsize',10, 'Color','w'); 
hcb.TickLabelInterpreter='latex';
hcb.FontSize = 11;
hcb.Color = 'w';
nexttile; hold on; xlabel('Distance from the center'), ylabel('Intensity (photons)');

plot(radius_ind, I_mean,'DisplayName', 'Average','LineWidth',2,'Color',[0 0.4470 0.7410]);
fill([radius_ind, fliplr(radius_ind)], [I_mean-I_std, fliplr(I_mean+I_std)], ...
    [0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.1);
set(gca,'FontSize',11,'Color','none','XColor','w', 'YColor','w');