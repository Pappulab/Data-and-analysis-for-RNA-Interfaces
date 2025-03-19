clear; close all;
%% v5 traj Tracking - requires results from v5 DiffuseAnalysis
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 2;
fovNum = 1;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
Tiffpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');

FigResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\fig_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\TrajTracks\');
if ~exist(FigResultsSavePath, 'dir')
   mkdir(FigResultsSavePath)
end

MatResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\mat_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\');


%% Load estimated results
load(strcat(MatResultsSavePath,'data',num2str(dataNum),'_FoV',num2str(fovNum),'_all.mat'),"molecule_numel_all","loc_list_sorted_all");

loc_list = loc_list_sorted_all;
molecule_list = molecule_numel_all;

%%%%%%%%%%% choose your mol
mol_choice = 2329;
mol_ind = 5;

disp("Lasting Frames:");
disp(molecule_list(molecule_list(:,1)==mol_choice,7));

FrameN = molecule_list(molecule_list(:,1)==mol_choice,7);

FirstFrame = max(molecule_list(molecule_list(:,1)==mol_choice,5)-1,1);
LastFrame = min(molecule_list(molecule_list(:,1)==mol_choice,6)+1,20000);

Frames = [FirstFrame:LastFrame];

loc_list_curr = loc_list(loc_list(:,1)>=FirstFrame&loc_list(:,1)<=LastFrame,:);

loc_mol = loc_list_curr(loc_list_curr(:,mol_ind)==mol_choice,:);

%% Load Raw tiff Image 
tiffList = dir(fullfile(Tiffpath, '*.tif'));

addpath(Tiffpath);
tiff_info = imfinfo(tiffList.name); % return tiff structure, one element per image
RawImg = imread(tiffList.name, Frames(1)) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : FrameN+2 % or set a fixed frame size!!!!
    temp_tiff = imread(tiffList.name, Frames(ii));
    RawImg = cat(3 , RawImg, temp_tiff);
end
RawImg = double(RawImg);

%%
clear tiff_info temp_tiff tiffList estList

%%
imgSz = size(RawImg,1:2);
frameSz = size(RawImg,3);
cropImg = RawImg(:,1:imgSz(2)/2,:);%1:imgSz(2)/2

loc_list_curr(:,2:3)=loc_list_curr(:,2:3)/58.5+(imgSz(1)+1)/2;
loc_mol(:,2:3)=loc_mol(:,2:3)/58.5+(imgSz(1)+1)/2;

rangeX = [max(floor(min(loc_mol(:,2)))-2, 1), min(ceil(max(loc_mol(:,2)))+2, imgSz(1))];
rangeY = [max(floor(min(loc_mol(:,3)))-5, 1), min(ceil(max(loc_mol(:,3)))+5, imgSz(1))];


filter_loc_ind = find(loc_list_curr(:,2)<rangeX(2) &  loc_list_curr(:,2)>rangeX(1) & ...
    loc_list_curr(:,3)<rangeY(2) & loc_list_curr(:,3)>rangeY(1));
loc_list_curr = loc_list_curr(filter_loc_ind,:);

%% Validation using the most distant trajectory
plotDiffuse_frames(loc_list_curr,mol_ind,mol_choice,cropImg,rangeX,rangeY,Frames,dataNum,fovNum,FigResultsSavePath);
close all
%%
function plotDiffuse_stack(loc_list,mol_select,mol_ind,cropImg,dataNum,savepath)
%%%%%%%%%%% This function is not changed yet 12.09pm 03202024
loc_select = loc_list(loc_list(:,mol_ind)==mol_select,:);
I_sz = size(cropImg,1);
frame_first = min(loc_select(:,1));
frame_last = max(loc_select(:,1));
loc_select(:,2:3)=loc_select(:,2:3)/58.5+(I_sz+1)/2;
figure('Position',[100,100,700,300]);


for frameNum = max(frame_first-1,1):min(frame_last+1,size(cropImg,3))
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; axis off; set(gca,'YDir','normal'); colorbar; hold on; caxis([95,145])
    
    line([I_sz+0.5,I_sz+0.5],[0.5,I_sz+0.5],'LineWidth',2,'color','w');
    if frameNum >= frame_first && frameNum<= frame_last
        idx = find(loc_select(:,1)==frameNum);
        if ~isempty(idx)
            SM_cur_loc = loc_select(idx,2:3);
            scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
            scatter(SM_cur_loc(1)+I_sz,SM_cur_loc(2),50,'rx','LineWidth',2);
            if frameNum==frame_first
                loc_temp = SM_cur_loc;
            elseif frameNum>=frame_first+1
                loc_temp = [loc_temp;SM_cur_loc]; 
            end
        end
               
        if frameNum>=frame_first+1
         line(loc_temp(:,1) ,loc_temp(:,2),'Color','k','LineWidth',2);  
        end
    end

    frame1 = getframe(gcf);
    im1 = frame2im(frame1);
    [A1,map1] = rgb2ind(im1,256);
    if frameNum == max(frame_first-1,1)
        imwrite(im1,[savepath,'\molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),...
            '.tif'],'tiff')
    else 
        imwrite(im1,[savepath,'\molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),...
            '.tif'],'tiff','WriteMode','append');
    end
    hold off
end
end
%%
function plotDiffuse_frames(loc_list_curr,mol_ind,mol_select,cropImg,rangeX,rangeY,imgFrames,dataNum,fovNum,savepath)
cropImg_new = cropImg(rangeY(1):rangeY(2),rangeX(1):rangeX(2),:);
loc_select = loc_list_curr;
loc_select(:,2) = loc_select(:,2)-rangeX(1)+1;
loc_select(:,3) = loc_select(:,3)-rangeY(1)+1;
mol_numel = unique(loc_select(:,mol_ind));
colorselect = [0 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330];
colorSz = size(colorselect,1);

% colormax = copper;
% colormax = colormax(end,:).^0.6;
% cMap  = [linspace(0,colormax(1),255)', linspace(0,colormax(2),255)', linspace(0,colormax(3),255)'];


for frameNum = 1:length(imgFrames)
    Fig = figure('Position',[600,490,323,375]);
    I_cur = (cropImg_new(:,:,frameNum)-100)*0.29;
    imagesc(I_cur(:,:)); axis image; axis off; 
    set(gca,'YDir','normal'); colorbar; colormap('copper'); hold on; caxis([8,25])

    for molNum = 1:length(mol_numel)    
        mol_curr = mol_numel(molNum);
        SM_loc = loc_select(loc_select(:,mol_ind)==mol_curr,:);
        
        SM_cur_loc = SM_loc(SM_loc(:,1)<=imgFrames(frameNum),:);
        sortrows(SM_cur_loc,1,'ascend');

        color_ind = mod(mol_curr-2,colorSz)+1;
%         color_ind(color_ind==0) = colorSz;
        color_inUse = colorselect(color_ind,:);
        if ~isempty(SM_cur_loc(SM_cur_loc(:,1)==imgFrames(frameNum),:))
            scatter(SM_cur_loc(end,2),SM_cur_loc(end,3),100,'x','MarkerEdgeColor',color_inUse,'LineWidth',2);
        end
        if size(SM_cur_loc,1)>1
            line(SM_cur_loc(:,2) ,SM_cur_loc(:,3),'Color',color_inUse,'LineWidth',1.5);  
            scatter(SM_cur_loc(1,2),SM_cur_loc(1,3),300,'.','MarkerEdgeColor',"w",'LineWidth',2);
        end
               

    end
    plot([size(cropImg_new,2)-500/58.5-0.5 size(cropImg_new,2)-0.5], [2, 2],'w','LineWidth', 8);

    exportgraphics(gcf,strcat(savepath,'data',num2str(dataNum),...
            '_FoV',num2str(fovNum),"_molecules_",num2str(mol_select),...
            "_frame",num2str(frameNum),'.pdf'),...
            'ContentType','vector','Resolution',600);
    hold off
end
end
