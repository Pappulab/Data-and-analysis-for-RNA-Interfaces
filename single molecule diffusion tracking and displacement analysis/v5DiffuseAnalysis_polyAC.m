clear; close all;
%% v4 diffuse analysis does not load raw images
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 8;
fovNum = 6;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');
ResultsSavePath =[fullpath,'\ResultsSave - 20240312'];
if ~exist(ResultsSavePath, 'dir')
   mkdir(ResultsSavePath)
end
%% Load estimated results
estList = dir(fullfile(fullpath, '*.mat'));
load(strcat(estList.folder,'\',estList.name));

%% Condition SMs
SM_est_save_all(:,1) = Angle_save(:,1);

SM_est_filtered = SM_est_save_all;
SM_est_filtered(:,2:3) = SM_est_filtered(:,2:3);
Angle_est_filtered = Angle_save;

%% filter outliers outside of condensate
r_thres_out= 1600;
r_thres_out= 2000;

Angle_est_filtered((SM_est_filtered(:,2)).^2+(SM_est_filtered(:,3)).^2>r_thres_out^2,:) = [];
SM_est_filtered((SM_est_filtered(:,2)).^2+(SM_est_filtered(:,3)).^2>r_thres_out^2,:) = [];

%% brightness filter
% add filtering here 
br_thres = 250;
Angle_est_filtered(SM_est_filtered(:,4)<br_thres,:) = [];
SM_est_filtered(SM_est_filtered(:,4)<br_thres,:) = [];

%% Angle filter
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

%% Check all loc data
figure('Color',[0,0,0],'Position',[460,290,375,435]); set(gca,'Color','w'); set(gcf,'color','none')
axis image; hold on;  title('all locs - eliminate outliers');
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'k','.'); 
% ax = gca;
% set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');
%% basic circle fit filter - for circle fit only!! Filter the centers
r_thres_base= 700;
x_center_base = -30; %150;
y_center_base = 0; %80;

SM_circFit = SM_est_filtered;
SM_circFit((SM_circFit(:,2)-x_center_base).^2+(SM_circFit(:,3)-y_center_base).^2<r_thres_base^2,:) = [];

%% Initial circular fit
% In fitting the curve, v2 updated using the nm scale coordinate (instead
% of pix) so the fitting circular is slightly different.
[R_init,XC_init,YC_init] = circfit(SM_circFit(:,2),SM_circFit(:,3));
plotcircfit(SM_circFit(:,2),SM_circFit(:,3));
hold on; scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'k','.'); 
%% Classify arcs
% indicate the classification bounds $arc1 = polyA-dilute, arc2 =  polyAC
% % 8 fov 1
% search_angle_arc1 = [-40,  140]; %$ dont set a bound right at 180 when you have to separate there
% search_angle_arc2 = [170, 180; -180, -60];
% % 8 fov 2
% search_angle_arc2 = [-150,  -30]; %$ dont set a bound right at 180 when you have to separate there
% search_angle_arc1 = [-20, 180; -180, -170];
% % 8 fov 3
% search_angle_arc1 = [-20,  165]; %$ dont set a bound right at 180 when you have to separate there
% search_angle_arc2 = [-170, -35];
% % 8 fov 4
% search_angle_arc1 = [150, 180; -180, -50]; %$ dont set a bound right at 180 when you have to separate there
% search_angle_arc2 = [-25,  115];
% % 8 fov 5
% search_angle_arc1 = [10,  160]; %$ dont set a bound right at 180 when you have to separate there
% search_angle_arc2 = [-175, -30];
% % 8 fov 6
search_angle_arc1 = [-60,  150]; %$ dont set a bound right at 180 when you have to separate there
search_angle_arc2 = [170, 180; -180, -80];

figure; hold on;
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
axis equal; title('all locs - with arc search bounds'); xlim([-2200,2200]),  ylim([-2200,2200])
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
for i = 1:size(search_angle_arc1,1)
    for j = 1:size(search_angle_arc1,2)
        if search_angle_arc1(i,j)==180 || search_angle_arc1(i,j)==-180
            continue
        else 
            line([XC_init,XC_init+(R_init+300)*cosd(search_angle_arc1(i,j))],...
                [YC_init,YC_init+(R_init+300)*sind(search_angle_arc1(i,j))],...
                'LineWidth',2,'color','r')
        end
    end
end
for i = 1:size(search_angle_arc2,1)
    for j = 1:size(search_angle_arc2,2)
        if search_angle_arc2(i,j)==180 || search_angle_arc2(i,j)==-180
            continue
        else 
            line([XC_init,XC_init+(R_init+300)*cosd(search_angle_arc2(i,j))],...
                [YC_init,YC_init+(R_init+300)*sind(search_angle_arc2(i,j))],...
                'LineWidth',2,'color','y')
        end
    end
end

% classify
SM_polar_angle = atan2d(SM_est_filtered(:,3)-YC_init,SM_est_filtered(:,2)-XC_init);
SM_arc1 = [];
SM_arc2 = [];
for i = 1:size(search_angle_arc1,1)
    ind = find(SM_polar_angle>search_angle_arc1(i,1) & SM_polar_angle<search_angle_arc1(i,2));
    SM_arc1 = [SM_arc1;SM_est_filtered(ind,:)];
end
for i = 1:size(search_angle_arc2,1)
    ind = find(SM_polar_angle>search_angle_arc2(i,1) & SM_polar_angle<search_angle_arc2(i,2));
    SM_arc2 = [SM_arc2;SM_est_filtered(ind,:)];
end
figure;
scatter(SM_arc1(:,2),SM_arc1(:,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
axis equal; title('arc 1 locs - poly A - dilute interface');  xlim([-2200,2200]),  ylim([-2200,2200])
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
figure;
scatter(SM_arc2(:,2),SM_arc2(:,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
axis equal; title('arc 2 locs - poly A/C interface');  xlim([-2200,2200]),  ylim([-2200,2200])
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
%% select arcs for next analysis (Filter fitting the new arc)
%%
% SM_est_filtered = SM_arc1;
% SM_circFit = SM_est_filtered;
% 
% % for arc 1 - poly A interface
% r_thres_base = 900;
% SM_circFit((SM_circFit(:,2)-XC_init).^2+(SM_circFit(:,3)-YC_init).^2<r_thres_base^2,:) = [];
% [R,XC,YC] = circfit(SM_circFit(:,2),SM_circFit(:,3));
% plotcircfit(SM_circFit(:,2),SM_circFit(:,3));
% hold on;
% scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),3,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
%% 
SM_est_filtered = SM_arc2;
SM_circFit = SM_est_filtered;
% for arc 2 - poly A/C interface
% LineVec1 = [0, -1]; LineVec2 = [-1,0]; LineVec3 = [-1,-1]; % data 8 FoV 1
% LineVec1 = [1, -1]; LineVec2 = [0,-1]; LineVec3 = [-1,-1]; % data 8 FoV 2
% LineVec1 = [0, -1]; LineVec2 = [1,-5]; LineVec3 = [-1,-1]; % data 8 FoV 3
% LineVec1 = [1, 0]; LineVec2 = [0,1]; LineVec3 = [1,1]; % data 8 FoV 4
% LineVec1 = [0, -1]; LineVec2 = [1,-5]; LineVec3 = [-1,-1]; % data 8 FoV 5
LineVec1 = [0, -1]; LineVec2 = [-1,0]; LineVec3 = [-1,-1]; % data 8 FoV 6

LineVec1 = LineVec1./norm(LineVec1); LineVec2 = LineVec2./norm(LineVec2); LineVec3 = LineVec3./norm(LineVec3);
Proj1 = sort(SM_circFit(:,2:3)*LineVec1','descend');
Proj2 = sort(SM_circFit(:,2:3)*LineVec2','descend'); 
Proj3 = sort(SM_circFit(:,2:3)*LineVec3','descend');
Point1 = mean(Proj1(1:round(0.1*length(Proj1))))*LineVec1;
Point2 = mean(Proj2(1:round(0.1*length(Proj2))))*LineVec2;
Point3 = mean(Proj3(1:round(0.1*length(Proj3))))*LineVec3;
ThreePoints = [Point1;Point2;Point3];
[R_tmp,XC_tmp,YC_tmp] = circfit(ThreePoints(:,1),ThreePoints(:,2));
clear LineVec1 LineVec2 LineVec3 Proj1 Proj2 Proj3 Points1 Points2 Points3

% r_thres_base = 900;

SM_circFit = SM_est_filtered;
SM_circFit((SM_circFit(:,2)-(XC_tmp)).^2+(SM_circFit(:,3)-(YC_tmp)).^2<(0.9*R_tmp)^2,:) = []; % for arc 2
[R,XC,YC] = circfit(SM_circFit(:,2),SM_circFit(:,3));
plotcircfit(SM_circFit(:,2),SM_circFit(:,3));
hold on;
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),3,'k','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
%% Filter the center region in this arc segmentation set up Because the center parts are meaningless!!!!!
% % or not, filter it after tracking
% % Use 0.36 R to filter
SM_est_filtered((SM_est_filtered(:,2)-(XC)).^2+(SM_est_filtered(:,3)-(YC)).^2<(0.5*R)^2,:) = [];
figure('Color',[0,0,0],'Position',[460,290,375,435]); set(gca,'Color','k'); set(gcf,'color','none')
axis image; hold on;  title('locs for diffusion analysis');
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'w','.'); 
%% Classify the loc into molecules
loc_list = SM_est_filtered(:,1:4);
loc_list(:,2:3) = loc_list(:,2:3);
frameN = max(loc_list(:,1));
threshold = 200;
check_frames = 2;

[loc_list, mol_ind] = classifyDiffuse(loc_list,threshold,check_frames,frameN);
%%
%%% Curent version 1: frame; 2-4: x, y, br; 5: mol ind; 6: traj length
%%% 7: mol interface angle
%%% 8: displacement; 
%%% 9: displacement azimuthal angle (average should be close to 0)
%%% 10-11: par/perp displacement
%%% Burst time of the diffusion analysis %% Adding distance
loc_list_sorted = sortrows(loc_list,mol_ind);
mol_range = 1:max(loc_list_sorted(:,mol_ind));
traj_length_frame = size(mol_range);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
dist_ind = size(loc_list_sorted,2);

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
for idx = 1:max(loc_list_sorted(:,mol_ind))
    traj_length_frame(idx) = max(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))-min(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,1))+1;
    if traj_length_frame(idx)>1
        loc_temp =  sortrows(loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,:),1);
        %loc_temp =  loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,:);
        % average jump (no angles)
        dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,dist_ind);
        
        for jj = 2:size(loc_temp,1)
            dist_temp = sqrt((loc_temp(jj-1,2)-loc_temp(jj,2))^2+(loc_temp(jj-1,3)-loc_temp(jj,3))^2);
            dist(jj) = dist(jj-1)+dist_temp;
           
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,dist_ind) = dist;
 
        total_dist = [total_dist;dist(end)];


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
            seq_jump_par(jj) = -[-vec_pos_norm(2),vec_pos_norm(1)]*vec_jump;
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_dist_ind) = seq_jump_dist;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_angle_ind) = seq_jump_angle;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_par_ind) = seq_jump_par;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,seq_jump_perp_ind) = seq_jump_perp;
    end
end
% jump histograms
%
jump_list = loc_list_sorted(:,seq_jump_dist_ind);
jump_list = jump_list(jump_list~=0);
jump_par_list = loc_list_sorted(:,seq_jump_par_ind);
jump_par_list = jump_par_list(jump_par_list~=0);
jump_perp_list = loc_list_sorted(:,seq_jump_perp_ind);
jump_perp_list = jump_perp_list(jump_perp_list~=0);
jump_angle_local = atan2d(jump_perp_list,jump_par_list); 
%% JUMP historgram filtering
% % (alternative jump historgram with filtering)
% RfilterRange = [666, 1000];
% mol_R = sqrt((loc_list_sorted(:,2)-XC).^2+(loc_list_sorted(:,3)-YC).^2);
% destinationFilter = find(mol_R>=RfilterRange(1)&mol_R<RfilterRange(2));
% jump_list = loc_list_sorted(destinationFilter,seq_jump_dist_ind);
% jump_list = jump_list(jump_list~=0);
% jump_par_list = loc_list_sorted(destinationFilter,seq_jump_par_ind);
% jump_par_list = jump_par_list(jump_par_list~=0);
% jump_perp_list = loc_list_sorted(destinationFilter,seq_jump_perp_ind);
% jump_perp_list = jump_perp_list(jump_perp_list~=0);
% jump_angle_local = atan2d(jump_perp_list,jump_par_list); 

%% Burst time and traj dist
burstTime = (traj_length_frame-1)*10;
figure('Position',[100,100,900,600]); subplot(211);
histogram(burstTime,'BinWidth',10.001);xlabel('Burst Time (ms)'), ylabel('count'); title(['Burst time among tracked molecules, threshold = ', num2str(threshold), 'nm']), xlim([0,250])
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
subplot(212); histogram(total_dist,100); xlabel('Distances (nm)'), ylabel('count'); title(['Total travel distance of molecules, threshold = ', num2str(threshold), 'nm']), xlim([0,1500])
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
%% Displacement and decomposed displacement
binsNum = 40; binsLimit = [-200,200];
figure('Position',[100,100,900,600]); 
subplot(211); histogram(jump_list,'Normalization','pdf','BinWidth',5);  
xlabel('Displacement per 10ms'), ylabel('pdf'); 
title(['Displacement per frame of molecules, threshold = ', num2str(threshold), 'nm'])
subtitle(['check ',num2str(check_frames),' frames']); xlim([0,400]);
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');set(ax.Subtitle,'Color','w');
subplot(212); hold on;
h_perp = histogram(jump_perp_list,binsNum,'BinLimits',binsLimit,'Normalization','pdf');
h_perp.DisplayName='perpendicular'; 
h_par = histogram(jump_par_list,binsNum,'BinLimits',binsLimit,'Normalization','pdf'); 
h_par.DisplayName='parallel'; 
xlabel('Displacement per 10ms'), ylabel('pdf');  ylim([0,0.01]), legend();
title('Displacement per frame of molecules'); xlim([-300,300]);
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
%% 1D gauss fit
% [pd_perp, pd_par, ci99_perp, ci99_par,figure1] = diffusion1Dgauss(jump_perp_list,jump_par_list,binsNum,binsLimit);
pd_perp = fitdist(jump_perp_list,'Normal');
pd_par = fitdist(jump_par_list,'Normal');
ci99_perp = paramci(pd_perp,'Alpha',.01);
ci99_par = paramci(pd_par,'Alpha',.01);
figure1 = figure('Position',[100,100,900,350]); 
h1 = histfit_YQdefined(jump_perp_list,"normal",binsLimit,binsNum);  hold on;
h2 = histfit_YQdefined(jump_par_list,"normal",binsLimit,binsNum); 
set(h1(1),'facecolor',"#D95319",'faceAlpha',0.7,'DisplayName','perpendicular'); set(h1(2),'color',"#D95319")
set(h2(1),'facecolor',"#4DBEEE",'faceAlpha',0.6,'DisplayName','parallel'); set(h2(2),'color', "#4DBEEE")
legend([h1(1),h2(1)]);

xlabel('Displacement per 10ms'), ylabel('counts');   
title('Displacement per frame of molecules'); xlim([-300,300]);

createtextbox_gaussfit_annotate(figure1,pd_perp,pd_par,ci99_perp,ci99_par);

ax = gca;
% set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');

% exportgraphics(gcf,strcat(ResultsSavePath,'\data_', num2str(dataNum),'FoV_', num2str(fovNum),'1DGauss_PolyA_C.pdf'),'Resolution',600);
%% Displacement angle in local coord
% figure('Position',[100,100,600,300]); 
% histogram(jump_angle_local,'Normalization','probability','BinWidth',6);  
% xlabel('Displacement angle in local coord'), ylabel('pdf');  ylim([0,0.05])
% title('Displacement angle per frame of molecules'); xlim([-180,180]);
%% Pol plot of the displacements and displacement angle
figure('Position',[100,100,600,300]); %sgtitle('Diffusion in local coordinates','Color','w')
subplot(121); polarhistogram((jump_angle_local)/180*pi, 36); title('Angles histogram');
% set(gca,'ThetaColor','w','RColor','k'); 
% ax = gca;
% set(ax,'Color','k','RColor','w','ThetaColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');
subplot(122);  
polarscatter((jump_angle_local)/180*pi,jump_list,4,'o','filled','MarkerFaceAlpha',0.3,'MarkerEdgeColor',"none"); 
% polarscatter((jump_angle_local)/180*pi,jump_list,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
title('Displacement');
% set(gca,'ThetaColor','w','RColor','k'); 
% ax = gca;
% set(ax,'Color','k','RColor','w','ThetaColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');

% exportgraphics(gcf,strcat(ResultsSavePath,'\data_', num2str(dataNum),'FoV_', num2str(fovNum),'polarplots_PolyA_C.pdf'),'Resolution',600);
%% 2D Gaussian fit
% [mu_x, mu_y, sigma_x, sigma_y, angle, fiterr, rr, fig1, fig2] = diffusion2Dgauss(jump_par_list,jump_perp_list);  

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
function plotDiffuse(loc_list,mol_select,mol_ind,cropImg,dataNum,threshold,check_frames,savepath)
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
        imwrite(im1,[savepath,'\threshold ', num2str(threshold),' checkFrames ', ...
            num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),...
            '.tif'],'tiff')
    else 
        imwrite(im1,[savepath,'\threshold ', num2str(threshold),' checkFrames ', ...
            num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),...
            '.tif'],'tiff','WriteMode','append');
    end
    hold off
end
end
%%
function plotTraj(loc_list,mol_select,mol_ind,dataNum,threshold,check_frames,path,colors,preNaming,axisRang)
colorSz = size(colors,1);
Fig1 = figure('Position',[475,114,740,600]); hold on; axis image
title(strcat("Trajtory - ",preNaming));
% subtitle(strcat('Molecules [', num2str(mol_select),']'));
xlabel('X position nm'); ylabel('Y position nm');
% xlim([-(threshold*1.5+50),(threshold*1.5+50)]); ylim([-(threshold*1.5+50),(threshold*1.5+50)]);
set(gca,'FontSize',14)

for i = 1:length(mol_select)
    mol_selectNum = mol_select(i);
    loc_select = loc_list(loc_list(:,mol_ind)==mol_selectNum,:);
    x_loc = loc_select(:,2);
    y_loc = loc_select(:,3);
    color_ind = mod(i,colorSz);
    color_ind(color_ind==0) = colorSz;
    colorSelect =  colors(color_ind,:);
    plot(x_loc, y_loc,'LineWidth',1.2,'Color',colorSelect); hold on;
    scatter(x_loc(1,:), y_loc(1,:),15,colorSelect(1:3),'filled');
end

if axisRang~=[]
    xlim([axisRang(1), axisRang(2)]);
    ylim([axisRang(3), axisRang(4)]);
else
    xlim([-1600, 1600]);
    ylim([-1600, 1600]);
end
exportgraphics(Fig1,strcat(path,'\threshold ', num2str(threshold),' checkFrames ', ...
            num2str(check_frames),' molecules',preNaming,'trajactory locs ', num2str(dataNum),'.jpg'),'Resolution',600);
end