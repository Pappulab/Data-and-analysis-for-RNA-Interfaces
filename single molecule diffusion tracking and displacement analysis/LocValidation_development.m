clear; close all;
dataNum = 48;
fullpath = pwd + "\data"+'\082423 YOPRO\'+num2str(dataNum);
%% Load Raw tiff Image and estimated results
tiffList = dir(fullfile(fullpath, '*.tif'));
estList = dir(fullfile(fullpath, '*.mat'));

addpath(fullpath);
tiff_info = imfinfo(tiffList.name); % return tiff structure, one element per image
RawImg = imread(tiffList.name, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(tiffList.name, ii);
    RawImg = cat(3 , RawImg, temp_tiff);
end
RawImg = double(RawImg);

load(estList.name);

clear tiff_info temp_tiff tiffList estList
%%
imgSz = size(RawImg,1:2);
frameSz = size(RawImg,3);
cropImg = RawImg(:,1:imgSz(2)/2,:);

%% Filtering
SM_est_save_all(:,1) = Angle_save(:,1);

SM_est_pix_filtered = SM_est_save_all;
SM_est_pix_filtered(:,2:3) = SM_est_pix_filtered(:,2:3)/58.5+45;
Angle_est_pix_filtered = Angle_save;
% YQ note: here exists a problem of pixel shifting, might caused by
% somewhere in RoSE-O estimation -> the center of the FoV is actually the
% (sz-1)/2


% add filtering here
br_thres=500;
Angle_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];
SM_est_pix_filtered(SM_est_pix_filtered(:,4)<br_thres,:) = [];

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
frameNum = 11;
% sample_frameSz = 200;
% for frameNum = 1:frameSz
SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)==frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3);
imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145])
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
%scatter(0+45,0+45,50,'rx','LineWidth',2);

% frame1 = getframe(gcf);
% im1 = frame2im(frame1);
% [A1,map1] = rgb2ind(im1,256);
% if frameNum == 1
%     imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),'.tif'],'tiff')
% else 
%     imwrite(im1,['First ',num2str(sample_frameSz), ' raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
% end
% end

%% Classify the loc into molecules
loc_list = SM_est_pix_filtered;
  loc_list(:,2:3) = (loc_list(:,2:3)-45)*58.5;
  threshold = 4*58.5;
  check_frames = 4;


loc_list_dim = size(loc_list);
mol_ind = loc_list_dim(2)+1;
% multiMol_ind = loc_list_dim(2)+2;
frameN = max(loc_list(:,1));
loc_list = [loc_list, zeros(loc_list_dim(1),1)];%,zeros(loc_list_dim(1),1)];
molNum = 0;

for frameInd = 1:frameN-1
    frame1 = frameInd;
    loc_f1 = loc_list(loc_list(:,1)==frame1,2:3);
    % giving zero idx molecules index
    mol_f1 = loc_list(loc_list(:,1)==frame1,mol_ind);
%     multiMol_1 = loc_list(loc_list(:,1)==frame1,multiMol_ind);
    for idx = 1:length(mol_f1)
        if mol_f1(idx)==0
            mol_f1(idx)=molNum+1;
            molNum = molNum+1;
        else
        end
    end
    loc_list(loc_list(:,1)==frame1,mol_ind) = mol_f1;
%     loc_list(loc_list(:,1)==frame1,multiMol_ind) = multiMol_1;
     
    for frameCur = 1:min(check_frames-1,frameN-frameInd)   
        frame_cur = frameInd+frameCur;
        loc_cur = loc_list(loc_list(:,1)==frame_cur,2:3);
        mol_cur = loc_list(loc_list(:,1)==frame_cur,mol_ind);
%         multiMol_cur = loc_list(loc_list(:,1)==frame_cur,multiMol_ind);
        
        [CurInd, f1Ind, d] = crossDistance(loc_f1,loc_cur,threshold);
        mol_cur(CurInd) = mol_f1(f1Ind);
%         multiMol_1(f1Ind) = 1;
%         multiMol_cur(CurInd) = 1;
        loc_list(loc_list(:,1)==frame_cur,mol_ind) = mol_cur;
%         loc_list(loc_list(:,1)==frame1,multiMol_ind) = multiMol_1;
%         loc_list(loc_list(:,1)==frame_cur,multiMol_ind) = multiMol_cur;
    end
end

if loc_list(loc_list(:,1)==frameN,mol_ind)==0
    loc_list(loc_list(:,1)==frameN,mol_ind)=molNum+1;
end

%% Burst time of the diffusion analysis %% Adding distance
loc_list_sorted = sortrows(loc_list,mol_ind);
mol_range = 1:max(loc_list_sorted(:,mol_ind));
traj_length = size(mol_range);
loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
dist_ind = size(loc_list_sorted,2);
loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
avg_dist_ind = size(loc_list_sorted,2);
total_dist=[];
total_avg_dist=[];
loc_list_sorted_allT = loc_list_sorted;
for idx = 1:max(loc_list_sorted(:,mol_ind))
    traj_length(idx) = max(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))-min(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))+1;
    if traj_length(idx)>1
        loc_temp = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,:);
        dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,dist_ind);
        avg_dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,avg_dist_ind);
        for jj = 2:size(loc_temp,1)
            dist_temp = sqrt((loc_temp(jj-1,2)-loc_temp(jj,2))^2+(loc_temp(jj-1,3)-loc_temp(jj,3))^2);
            dist(jj) = dist(jj-1)+dist_temp;
            avg_dist(jj) = dist(jj)/(loc_temp(jj,1)-loc_temp(1,1));
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,dist_ind) = dist;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,avg_dist_ind) = avg_dist;
        total_dist = [total_dist;dist(end)];
        total_avg_dist = [total_avg_dist;avg_dist(end)];
    end
end
%%
figure('Position',[100,100,600,400]); histogram(traj_length(traj_length>1)*10); xlabel('Burst Time (ms)'), ylabel('counts'); title('Burst time among tracked molecules')
figure('Position',[100,100,600,400]); histogram(total_dist); xlabel('Distances (nm)'), ylabel('counts'); title('Total travel distance of molecules')
figure('Position',[100,100,600,400]); histogram(total_avg_dist);  xlabel('Travel speed (nm/10ms)'), ylabel('counts'); title('Travel distances per frame of molecules')
%% Validation using the longest trajectory
[l_length,l_idx] = max(traj_length);

loc_longest = loc_list_sorted(loc_list_sorted(:,mol_ind)==l_idx,:);
frame_first = min(loc_longest(:,1));
frame_last = max(loc_longest(:,1));
loc_longest(:,2:3)=loc_longest(:,2:3)/58.5+45;

for frameNum = frame_first:frame_last
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145])
    idx = find(loc_longest(:,1)==frameNum);
    if ~isempty(idx)
        SM_cur_loc = loc_longest(idx,2:3);
        scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
        if frameNum==frame_first
            loc_temp = SM_cur_loc;
        elseif frameNum>=frame_first+1
            loc_temp = [loc_temp;SM_cur_loc]; 
        end

%scatter(0+45,0+45,50,'rx','LineWidth',2);
    else
    end
if frameNum>=frame_first+1
 line(loc_temp(:,1) ,loc_temp(:,2),'Color','#D95319','LineWidth',2);  
end
frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
if frameNum == frame_first
    imwrite(im1,['Longest trajactory molecule ', num2str(l_idx), ' raw and locs ', num2str(dataNum),'.tif'],'tiff')
else 
    imwrite(im1,['Longest trajactory molecule ', num2str(l_idx), ' raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
end
end
%% Validation using the most distant trajectory
[f_length,f_idx] = max(loc_list_sorted(:,avg_dist_ind));
f_mol_ind = loc_list_sorted(f_idx,mol_ind);
loc_dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==f_mol_ind,:);
frame_first = min(loc_dist(:,1));
frame_last = max(loc_dist(:,1));
loc_dist(:,2:3)=loc_dist(:,2:3)/58.5+45;

for frameNum = frame_first:frame_last
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145])
    idx = find(loc_dist(:,1)==frameNum);
    if ~isempty(idx)
        SM_cur_loc = loc_dist(idx,2:3);
        scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
        if frameNum==frame_first
            loc_temp = SM_cur_loc;
        elseif frameNum>=frame_first+1
            loc_temp = [loc_temp;SM_cur_loc]; 
        end

%scatter(0+45,0+45,50,'rx','LineWidth',2);
    else
    end
if frameNum>=frame_first+1
 line(loc_temp(:,1) ,loc_temp(:,2),'Color','#D95319','LineWidth',2);  
end
frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
if frameNum == frame_first
    imwrite(im1,['Furtherest trajactory molecule ', num2str(f_mol_ind), ' raw and locs ', num2str(dataNum),'.tif'],'tiff')
else 
    imwrite(im1,['Furtherest trajactory molecule ', num2str(f_mol_ind), ' raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
end
end


%% Validation using the most distant trajectory
mol_select = 175; %66
loc_select = loc_list_sorted(loc_list_sorted(:,mol_ind)==mol_select,:);
frame_first = min(loc_select(:,1));
frame_last = max(loc_select(:,1));
loc_select(:,2:3)=loc_select(:,2:3)/58.5+45;

for frameNum = frame_first:frame_last
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145])
    idx = find(loc_select(:,1)==frameNum);
    if ~isempty(idx)
        SM_cur_loc = loc_select(idx,2:3);
        scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
        if frameNum==frame_first
            loc_temp = SM_cur_loc;
        elseif frameNum>=frame_first+1
            loc_temp = [loc_temp;SM_cur_loc]; 
        end

%scatter(0+45,0+45,50,'rx','LineWidth',2);
    else
    end
if frameNum>=frame_first+1
 line(loc_temp(:,1) ,loc_temp(:,2),'Color','#D95319','LineWidth',2);  
end
frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
if frameNum == frame_first
    imwrite(im1,['molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff')
else 
    imwrite(im1,['molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
end
end
%% 
function [Vind, Uind, d] = crossDistance(U,V,threshold)
% U, V must be i-by-2 and j-by-2 x-y coordinated locations
if ~isempty(U) && ~isempty(V)
    Ux = U(:,1); Uy = U(:,2);
    Vx = V(:,1); Vy = V(:,2);
    [UX,VX] = meshgrid(Ux,Vx);
    [UY,VY] = meshgrid(Uy,Vy);
    dMat = sqrt((UX-VX).^2+(UY-VY).^2);
    [Vind, Uind, d] = find(dMat<=threshold);
else
    Vind = [];
    Uind = [];
    d = [];
end
end