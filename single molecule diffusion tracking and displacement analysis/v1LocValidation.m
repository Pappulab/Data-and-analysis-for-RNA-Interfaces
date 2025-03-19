clear; close all;
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 12;
fovNum = 1;
resultFolder = '\20240207 test new Nanoscope';
dataFolder = '\data';
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');

%% Load Raw tiff Image
tiffList = dir(fullfile(fullpath, '*.tif'));

addpath(fullpath);
tiff_info = imfinfo(tiffList.name); % return tiff structure, one element per image
RawImg = imread(tiffList.name, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(tiffList.name, ii);
    RawImg = cat(3 , RawImg, temp_tiff);
end
RawImg = double(RawImg);

%% Load estimated results
estList = dir(fullfile(fullpath, '*.mat'));
load(strcat(estList.folder,'\',estList.name));
%%
clear tiff_info temp_tiff tiffList estList
%%
imgSz = size(RawImg,1:2);
frameSz = size(RawImg,3);
cropImg = RawImg(:,:,:);

%% Filtering
SM_est_save_all(:,1) = Angle_save(:,1); % for 2D
% SM_est_save_all = [Angle_save(:,1),SM_est_save_all]; % for 3D
%%
SM_est_pix_filtered = SM_est_save_all;
SM_est_pix_filtered(:,2:3) = SM_est_pix_filtered(:,2:3)/58.5+(imgSz(1)+1)/2;
Angle_est_pix_filtered = Angle_save;
% YQ note: here exists a problem of pixel shifting, might caused by
% somewhere in RoSE-O estimation -> the center of the FoV is actually the
% (sz-1)/2
%%

% add filtering here
br_thres=400;
% For 2D
Angle_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];
SM_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];
% For 3D
% Angle_est_pix_filtered(SM_est_pix_filtered(:,5)<br_thres,:) = [];
% SM_est_pix_filtered(SM_est_pix_filtered(:,5)<br_thres,:) = [];

% Angle_est_pix_filtered(SM_est_pix_filtered(:,8)<-0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,8)<-0.75,:) =[];
% Angle_est_pix_filtered(SM_est_pix_filtered(:,8)>0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,8)>0.75,:) =[];
% Angle_est_pix_filtered(SM_est_pix_filtered(:,9)<-0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,9)<-0.75,:) =[];
% Angle_est_pix_filtered(SM_est_pix_filtered(:,9)>0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,9)>0.75,:) =[];
% Angle_est_pix_filtered(SM_est_pix_filtered(:,10)<-0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,10)<-0.75,:) =[];
% Angle_est_pix_filtered(SM_est_pix_filtered(:,10)>0.75,:) =[];
% SM_est_pix_filtered(SM_est_pix_filtered(:,10)>0.75,:) =[];


%%
frameNum = 1;
% sample_frameSz = 200;
% for frameNum = 1:frameSz
SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)==frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3); figure('Position',[100,100,700,300]);
imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145]); colormap("copper");
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
scatter(SM_cur_loc(:,1)+imgSz(1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
% scatter(1,1,50,'rx','LineWidth',2)
% scatter(imgSz(1),1,50,'rx','LineWidth',2)
line([imgSz(1)+0.5,imgSz(1)+0.5],[0,imgSz(1)],'LineWidth',2,'color','w');
line([imgSz(1)*2-14,imgSz(1)*2-14+8.5470],[imgSz(1)-4,imgSz(1)-4],'LineWidth',5,'color','w');
%scatter(0+45,0+45,50,'rx','LineWidth',2);
% set(gca,'YDir','normal')
% frame1 = getframe(gcf);
% im1 = frame2im(frame1);
% [A1,map1] = rgb2ind(im1,256);
% if frameNum == 1
%     imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),'.tif'],'tiff')
% else 
%     imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
% end
% end
% set(gcf, 'InvertHardcopy', 'on');
% whitebg('k');
%%
%%
sample_frameSz = 100;
for frameNum = 1:sample_frameSz
SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)==frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3); figure('Position',[100,100,700,300]);
imagesc(cropImg(:,:,frameNum)); axis image; axis off; set(gca,'YDir','normal'); colorbar; hold on; caxis([95,145])
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
scatter(SM_cur_loc(:,1)+imgSz(1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
line([imgSz(1)+0.5,imgSz(1)+0.5],[0,imgSz(1)],'LineWidth',2,'color','w');
line([imgSz(1)*2-14,imgSz(1)*2-14+8.5470],[imgSz(1)-4,imgSz(1)-4],'LineWidth',5,'color','w');
%scatter(0+45,0+45,50,'rx','LineWidth',2);

frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
if frameNum == 1
    imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),' total loc ', num2str(size(SM_est_pix_filtered,1)),'.tif'],'tiff')
else 
    imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),' total loc ', num2str(size(SM_est_pix_filtered,1)),'.tif'],'tiff','WriteMode','append');
end
end
% set(gcf, 'InvertHardcopy', 'on');
% whitebg('k');