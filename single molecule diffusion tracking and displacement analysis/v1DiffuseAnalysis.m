clear; close all;
%%
dataNum = 48;
resultFolder = '\20231019 Results\';
dataFolder = '\101923 YOPRO\';
fullpath = strcat(pwd, resultFolder, "data",dataFolder,num2str(dataNum));
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
% clearvars -except RawImg SM_est_save_all Angle_save dataNum
%%
imgSz = size(RawImg,1:2);
frameSz = size(RawImg,3);
cropImg = RawImg(:,:,:);%1:imgSz(2)/2

%% Filtering
SM_est_save_all(:,1) = Angle_save(:,1);

SM_est_pix_filtered = SM_est_save_all;
SM_est_pix_filtered(:,2:3) = SM_est_pix_filtered(:,2:3)/58.5+45;
Angle_est_pix_filtered = Angle_save;
% YQ note: here exists a problem of pixel shifting, might caused by
% somewhere in RoSE-O estimation -> the center of the FoV is actually the
% (sz-1)/2

%%
% add filtering here
br_thres=000;
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

%% Circule fit
[R,XC,YC] = circfit(SM_est_pix_filtered(:,2),SM_est_pix_filtered(:,3));
plotcircfit(SM_est_pix_filtered(:,2),SM_est_pix_filtered(:,3));
%% check frame
frameNum = 19;
% sample_frameSz = 200;
% for frameNum = 1:frameSz
SM_cur_loc = [];
idx = find(SM_est_pix_filtered(:,1)==frameNum);
SM_cur_loc = SM_est_pix_filtered(idx,2:3);
imagesc(cropImg(:,:,frameNum)); axis image; colorbar; hold on; caxis([95,145])
scatter(SM_cur_loc(:,1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
scatter(SM_cur_loc(:,1)+imgSz(1),SM_cur_loc(:,2),50,'rx','LineWidth',2);
line([imgSz(1)+0.5,imgSz(1)+0.5],[0,imgSz(1)],'LineWidth',2,'color','w');

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
loc_list = SM_est_pix_filtered(:,1:4);
loc_list(:,2:3) = (loc_list(:,2:3)-45)*58.5;
frameN = max(loc_list(:,1));
threshold = 200;
check_frames = 2;

[loc_list, mol_ind] = classifyDiffuse(loc_list,threshold,check_frames,frameN);

%%% Burst time of the diffusion analysis %% Adding distance
%%% Curent version 1: frame; 2-4: x, y, br; 5: mol ind; 6: traj length
%%% 7: average displace per traj; 8: mol interface angle
%%% 9: displacement; 10: displacement azimuthal angle
%%% 11-12: par/perp displacement
loc_list_sorted = sortrows(loc_list,mol_ind);
mol_range = 1:max(loc_list_sorted(:,mol_ind));
traj_length_frame = size(mol_range);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
dist_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
avg_dist_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
mol_angle_ind = size(loc_list_sorted,2);
loc_list_sorted(:,mol_angle_ind) = atan2d(loc_list_sorted(:,3)-YC,loc_list_sorted(:,2)-XC);
%(might work for evaluation but not used here or later yet)

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_dist_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_angle_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_par_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_perp_ind = size(loc_list_sorted,2);


total_dist=[];
total_avg_dist=[];
for idx = 1:max(loc_list_sorted(:,mol_ind))
    traj_length_frame(idx) = max(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))-min(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))+1;
    if traj_length_frame(idx)>1
        loc_temp = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,:);
        % average jump (no angles)
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

        % sequential jumps (with angle)
        seq_jump_dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_dist_ind);
        seq_jump_angle = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_angle_ind);
        seq_jump_par = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_par_ind);
        seq_jump_perp = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_perp_ind);
        for jj = 2:size(loc_temp,1)
            vec_pos = [loc_temp(jj,2)-XC,loc_temp(jj,3)-YC]';
            vec_pos_norm = vec_pos/sqrt(sum(vec_pos.^2));
            vec_jump = [loc_temp(jj,2)-loc_temp(jj-1,2),loc_temp(jj,3)-loc_temp(jj-1,3)]';
            jump_temp = sqrt((loc_temp(jj-1,2)-loc_temp(jj,2))^2+(loc_temp(jj-1,3)-loc_temp(jj,3))^2);
            seq_jump_dist(jj) = jump_temp/(loc_temp(jj,1)-loc_temp(jj-1,1));
            seq_jump_angle(jj) = atan2d(vec_jump(2),vec_jump(1));
            seq_jump_perp(jj) = vec_pos_norm'*vec_jump;
            if sign(vec_pos_norm(1)*vec_pos_norm(2))== 1
                seq_jump_par(jj) = -[-vec_pos_norm(2),vec_pos_norm(1)]*vec_jump;
            else
                seq_jump_par(jj) = [-vec_pos_norm(2),vec_pos_norm(1)]*vec_jump;
            end 
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_dist_ind) = seq_jump_dist;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_angle_ind) = seq_jump_angle;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_par_ind) = seq_jump_par;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_perp_ind) = seq_jump_perp;
    end
end
total_mol_num_not1fr = length(total_avg_dist);% molecule number that lasts more than 1 frame
%%
figure('Position',[100,100,600,400]); 
histogram(traj_length_frame(traj_length_frame>1)*10,40); xlabel('Burst Time (ms)'), ylabel('count'); title(['Burst time among tracked molecules, threshold = ', num2str(threshold), 'nm'])
figure('Position',[100,100,600,400]); histogram(total_dist,40); xlabel('Distances (nm)'), ylabel('count'); title(['Total travel distance of molecules, threshold = ', num2str(threshold), 'nm'])

%%
figure('Position',[100,100,600,400]); histogram(total_avg_dist,'Normalization','probability','BinWidth',5);  
xlabel('Travel speed (nm/10ms)'), ylabel('pdf'); 
title(['Travel distances per frame of molecules, threshold = ', num2str(threshold), 'nm'])
subtitle(['check ',num2str(check_frames),' frames'])
xlim([0,400]);
% ,'BinLimits',[0,100]


%% jump histograms
%%
jump_list = loc_list_sorted(:,seq_jump_dist_ind);
jump_list = jump_list(jump_list~=0);
jump_par_list = loc_list_sorted(:,seq_jump_par_ind);
jump_par_list = jump_par_list(jump_par_list~=0);
jump_perp_list = loc_list_sorted(:,seq_jump_perp_ind);
jump_perp_list = jump_perp_list(jump_perp_list~=0);

%%
figure('Position',[100,100,600,400]); 
histogram(jump_list,'Normalization','probability','BinWidth',5);  
xlabel('Displacement per 10ms'), ylabel('pdf'); 
title(['Displacement per frame of molecules, threshold = ', num2str(threshold), 'nm'])
subtitle(['check ',num2str(check_frames),' frames']); xlim([0,400]);
%%
binsNum = 40; binsLimit = [-200,200];
figure('Position',[100,100,900,600]); 
subplot(211); histogram(jump_list,'Normalization','probability','BinWidth',5);  
xlabel('Displacement per 10ms'), ylabel('pdf'); 
title(['Displacement per frame of molecules, threshold = ', num2str(threshold), 'nm'])
subtitle(['check ',num2str(check_frames),' frames']); xlim([0,400]);
subplot(212); hold on;
h_perp = histogram(jump_perp_list,binsNum,'BinLimits',binsLimit,'Normalization','probability');
h_perp.DisplayName='perpendicular'; 
h_par = histogram(jump_par_list,binsNum,'BinLimits',binsLimit,'Normalization','probability'); 
h_par.DisplayName='parallel'; 
xlabel('Displacement per 10ms'), ylabel('pdf');  ylim([0,0.1]), legend();
title('Displacement per frame of molecules'); xlim([-300,300]);
%% gauss fit
pd_perp = fitdist(jump_perp_list,'Normal');
pd_par = fitdist(jump_par_list,'Normal');
ci99_perp = paramci(pd_perp,'Alpha',.01);
ci99_par = paramci(pd_par,'Alpha',.01);
gass_plot = @(d,sigma,mu) 1/sqrt(2*pi*sigma^2)*exp(-1/2*((d-mu)^2/sigma^2));
figure1 = figure('Position',[100,100,900,350]); 
h1 = histfit_YQdefined(jump_perp_list,"normal",binsLimit,binsNum);  hold on;
h2 = histfit_YQdefined(jump_par_list,"normal",binsLimit,binsNum); 
set(h1(1),'facecolor',"#D95319",'faceAlpha',0.7,'DisplayName','perpendicular'); set(h1(2),'color',"#D95319")
set(h2(1),'facecolor',"#4DBEEE",'faceAlpha',0.6,'DisplayName','parallel'); set(h2(2),'color', "#4DBEEE")
legend([h1(1),h2(1)]);

xlabel('Displacement per 10ms'), ylabel('counts');   
title('Displacement per frame of molecules'); xlim([-300,300]);

createtextbox_gaussfit_annotate(figure1,pd_perp,pd_par,ci99_perp,ci99_par);
%%
jump_angle_local = atan2d(jump_perp_list,jump_par_list); 
figure('Position',[100,100,600,300]); 
histogram(jump_angle_local,'Normalization','probability','BinWidth',6);  
xlabel('Displacement angle in local coord'), ylabel('pdf');  ylim([0,0.05])
title('Displacement angle per frame of molecules'); xlim([-180,180]);
%%
figure('Position',[100,100,300,300]); polarscatter((jump_angle_local+180)/180*pi,jump_list,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); set(gca,'ThetaColor','w','RColor','k');
%% sliding window average
stepSz = 2;
pol_ind = -180:stepSz:180;
avg_jump_direction = zeros(size(pol_ind));
for i=1:length(pol_ind)
    jump_temp_ind = find(jump_angle_local>=pol_ind(i)-stepSz/2 & jump_angle_local<pol_ind(i)+stepSz/2);
    avg_jump_direction(i) = mean(jump_list(jump_temp_ind));
end
figure('Position',[100,100,300,300]);
polarplot((pol_ind+180)/180*pi,avg_jump_direction);set(gca,'ThetaColor','w','RColor','k');

%% 2D Gaussian fit
%%
figure('Position',[100,100,500,500]); scatter(jump_par_list,jump_perp_list,'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.5); axis equal;
    
%% Validation using the most distant trajectory
mol_select = 922; %66
plotDiffuse(loc_list,mol_select,mol_ind,cropImg,dataNum,threshold,check_frames)

%% Validation using the longest trajectory
[l_length,l_idx] = max(traj_length_frame);
plotDiffuse(loc_list,l_idx,mol_ind,cropImg,'longest',threshold,check_frames)

%% Validation using the most distant trajectory
[f_length,f_idx] = max(loc_list_sorted(:,dist_ind));
f_mol_ind = loc_list_sorted(f_idx,mol_ind);
plotDiffuse(loc_list,f_mol_ind,mol_ind,cropImg,'furthest',threshold,check_frames)

%%

%% 
function [Vind, Uind, d] = crossDistance(U,V,threshold)
% U, V must be i-by-2 and j-by-2 x-y coordinated locations
if ~isempty(U) && ~isempty(V)
    Ux = U(:,1); Uy = U(:,2);
    Vx = V(:,1); Vy = V(:,2);
    [UX,VX] = meshgrid(Ux,Vx);
    [UY,VY] = meshgrid(Uy,Vy);
    dMat = sqrt((UX-VX).^2+(UY-VY).^2);
    [Vind, Uind, d] = find(min(dMat,[],"all")&dMat<=threshold);
else
    Vind = [];
    Uind = [];
    d = [];
end
end

function [Vind, d] = crossDistance1D(U,V,threshold)
% U, V must be i-by-2 and j-by-2 x-y coordinated locations
if ~isempty(U) && ~isempty(V)
    Ux = U(:,1); Uy = U(:,2);
    Vx = V(:,1); Vy = V(:,2);
    [UX,VX] = meshgrid(Ux,Vx);
    [UY,VY] = meshgrid(Uy,Vy);
    dMat = sqrt((UX-VX).^2+(UY-VY).^2);
    [d,Vind] = min(dMat,[],'all');
    if d>threshold
        Vind = [];
        d = [];
    end
else
    Vind = [];
    d = [];
end
end
%%
function [loc_list,mol_ind] = classifyDiffuse(loc_list,threshold,check_frames,frameN)
loc_list_dim = size(loc_list);
mol_ind = loc_list_dim(2)+1;
% multiMol_ind = loc_list_dim(2)+2;
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
        
        for iii = 1:size(loc_f1,1)
            if ismember(mol_f1(iii),mol_cur)
                continue
            else
            [CurInd, d] = crossDistance1D(loc_f1(iii,:),loc_cur,threshold);
            mol_cur(CurInd) = mol_f1(iii);
            if ~isempty(CurInd)
                loc_f1(iii,:) = (loc_f1(iii,:)*frameCur+loc_cur(CurInd,:))/(frameCur+1);
            end
            end
        end

        % [CurInd, f1Ind, d] = crossDistance(loc_f1,loc_cur,threshold);
        % mol_cur(CurInd) = mol_f1(f1Ind);
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
end

%%
function plotDiffuse(loc_list,mol_select,mol_ind,cropImg,dataNum,threshold,check_frames)
loc_select = loc_list(loc_list(:,mol_ind)==mol_select,:);
frame_first = min(loc_select(:,1));
frame_last = max(loc_select(:,1));
loc_select(:,2:3)=loc_select(:,2:3)/58.5+45;
figure('Position',[100,100,700,300]);
I_sz = size(cropImg,1);

for frameNum = frame_first:frame_last
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; axis off; colorbar; hold on; caxis([95,145])
    idx = find(loc_select(:,1)==frameNum);
    if ~isempty(idx)
        SM_cur_loc = loc_select(idx,2:3);
        scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
        scatter(SM_cur_loc(1)+I_sz,SM_cur_loc(2),50,'rx','LineWidth',2);
        line([I_sz+0.5,I_sz+0.5],[0,I_sz],'LineWidth',2,'color','w');
        if frameNum==frame_first
            loc_temp = SM_cur_loc;
        elseif frameNum>=frame_first+1
            loc_temp = [loc_temp;SM_cur_loc]; 
        end

%scatter(0+45,0+45,50,'rx','LineWidth',2);
    else
    end
if frameNum>=frame_first+1
 line(loc_temp(:,1) ,loc_temp(:,2),'Color','k','LineWidth',2);  
end
frame1 = getframe(gcf);
im1 = frame2im(frame1);
[A1,map1] = rgb2ind(im1,256);
if frameNum == frame_first
    imwrite(im1,['threshold ', num2str(threshold),' checkFrames ', num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff')
else 
    imwrite(im1,['threshold ', num2str(threshold),' checkFrames ', num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
end
end
end