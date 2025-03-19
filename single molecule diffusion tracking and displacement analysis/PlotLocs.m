 clear; close all;
%% v5 diffuse analysis does not load raw images
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 2;
fovNum = 1;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');

Nimg = 1000;

%%
FigResultsSavePath = strcat("C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures\ProjectReview\fig_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\Loc_raw\');
if ~exist(FigResultsSavePath, 'dir')
   mkdir(FigResultsSavePath)
end

%% Load estimated results
estList = dir(fullfile(fullpath, '*.mat'));
load(strcat(estList.folder,'\',estList.name));

clear estList
%%
tiffList = dir(fullfile(fullpath, '*.tif'));

SMLMName = tiffList.name; % return tiff structure, one element per image
count =0;
for jj = 1:Nimg
        count = count+1;
        SMLMR = Tiff(SMLMName,'r');
        setDirectory(SMLMR,jj);
        SM_img(:,:,count) = (double(SMLMR.read)-100)*0.29;
end
%%
%%
imgSz = size(SM_img,1:2);
frameSz = size(SM_img,3);
cropImg = SM_img(:,:,:);

%% Condition SMs
SM_est_save_all(:,1) = Angle_save(:,1);

SM_est_filtered = SM_est_save_all;
SM_est_filtered(:,2:3) = SM_est_filtered(:,2:3);
Angle_est_filtered = Angle_save;
%%
SM_est_pix_filtered = SM_est_filtered;
SM_est_pix_filtered(:,2:3) = SM_est_pix_filtered(:,2:3)/58.5+(imgSz(1)+1)/2;
Angle_est_pix_filtered = Angle_save;

%%

% add filtering here
br_thres=400;
% For 2D
Angle_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];
SM_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];

%%
% frameNum = 100;
% % sample_frameSz = 200;
for frameNum = 1:frameSz
SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)==frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3); 


figure('Position',[439.8571428571428,299.2857142857143,1000,335.4285714285714]); 
subplot(121);

imagesc(cropImg(:,1:imgSz(1),frameNum)); axis image; hcb = colorbar; hold on; caxis([8,20]);  colormap("copper");
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2),400,[199 160 151]/255,'LineWidth',2,'MarkerFaceAlpha',0.9); 
line([imgSz(1)-25,imgSz(1)-25+8.*2],[imgSz(1)-4,imgSz(1)-4],'LineWidth',5,'color','w');

set(hcb.Title,'String','photons')
axis off;


SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)<=frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3); 

subplot(122);
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2), 25,[199 160 151]/255,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0); 

axis image;
xlim([0, imgSz(1)]); ylim([0, imgSz(1)])
set(gca,'Color','k');
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca, 'TickLength',[0 0])



set(gca,'YDir','normal')
frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
    if frameNum == 1
        imwrite(im1,strcat(FigResultsSavePath,'localization_raw.tif'),'tiff','Compression', 'none');
    else 
        imwrite(im1,strcat(FigResultsSavePath,'localization_raw.tif'),'tiff','Compression', 'none','WriteMode','append');
    end
close all
end
set(gcf, 'InvertHardcopy', 'on');
