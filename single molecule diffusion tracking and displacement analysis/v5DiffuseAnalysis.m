clear; close all;
%% v5 diffuse analysis does not load raw images
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 8;
fovNum = 1;
for i = 1:6
    fovNum=i;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');
% ResultsSavePath =[fullpath,'\ResultsSave - 20240312'];
% if ~exist(ResultsSavePath, 'dir')
%    mkdir(ResultsSavePath)
% end
%%
MatResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\mat_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\');
if ~exist(MatResultsSavePath, 'dir')
   mkdir(MatResultsSavePath)
end

FigResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\fig_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\');
if ~exist(FigResultsSavePath, 'dir')
   mkdir(FigResultsSavePath)
end

%% Set fonts and interpreter

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')
%% Load estimated results
estList = dir(fullfile(fullpath, '*.mat'));
load(strcat(estList.folder,'\',estList.name));

clear estList
%% Condition SMs
SM_est_save_all(:,1) = Angle_save(:,1);

SM_est_filtered = SM_est_save_all;
SM_est_filtered(:,2:3) = SM_est_filtered(:,2:3);
Angle_est_filtered = Angle_save;

Fig = figure('Position',[600,490,323,375]);  set(gca,'Color','k');
axis image; hold on;  title('all locs');
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'w','.'); 

% ax = gca;
% set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');
%
imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_All_locs');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%% filter outliers outside of condensate
r_thres_outlier= 1800;
SM_est_filtered((SM_est_filtered(:,2)).^2+(SM_est_filtered(:,3)).^2>r_thres_outlier^2,:) = [];

[~,X_tmp,Y_tmp] = circfit(SM_est_filtered(:,2),SM_est_filtered(:,3));

R_locs_tmp = sqrt((SM_est_filtered(:,2)-X_tmp).^2+(SM_est_filtered(:,3)-Y_tmp).^2);
R_locs_tmp_sorted = sort(R_locs_tmp,'descend');
R_out_threshold = 1.2*std(R_locs_tmp_sorted) ...
    + mean(R_locs_tmp_sorted(1:round(0.15*length(R_locs_tmp_sorted))));


Angle_est_filtered((SM_est_filtered(:,2)-X_tmp).^2+(SM_est_filtered(:,3)-Y_tmp).^2>R_out_threshold^2,:) = [];
SM_est_filtered((SM_est_filtered(:,2)-X_tmp).^2+(SM_est_filtered(:,3)-Y_tmp).^2>R_out_threshold^2,:) = [];

clear R_locs_tmp R_locs_tmp_sorted 
%% Check all loc data
Fig = figure('Position',[460,290,323,375]);  set(gca,'Color','k');
axis image; hold on;  title('all localizations - eliminate outliers');
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'w','.'); 
xlim([X_tmp-1.2*R_out_threshold,X_tmp+1.2*R_out_threshold]); ylim([Y_tmp-1.2*R_out_threshold,Y_tmp+1.2*R_out_threshold]);
hold on,
set(gca,'XTickLabel',[],'YTickLabel',[]);
plot([min(X_tmp+0.72*R_out_threshold,X_tmp+1.2*R_out_threshold-600),min(X_tmp+0.72*R_out_threshold,X_tmp+1.2*R_out_threshold-600)+500], [Y_tmp-1.1*R_out_threshold+100, Y_tmp-1.1*R_out_threshold+100],'w','LineWidth', 4);
% text(min(X_tmp+0.72*R_out_threshold,X_tmp+1.2*R_out_threshold-600)+250,Y_tmp-1.1*R_out_threshold+260,'500 nm','Color','w','FontSize',9,'FontWeight','bold','Interpreter','latex','HorizontalAlignment', 'center');
% ax = gca;
% set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');
imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_Filtered_all_locs');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

clear X_tmp Y_tmp  R_out_threshold

%% The first circle fit 
[~,X_locs,Y_locs] = circfit(SM_est_filtered(:,2),SM_est_filtered(:,3));

R_locs_tmp = sqrt((SM_est_filtered(:,2)-X_locs).^2+(SM_est_filtered(:,3)-Y_locs).^2);
R_angle_tmp = atan2d(SM_est_filtered(:,3)-Y_locs,SM_est_filtered(:,2)-X_locs);
[R_locs_tmp,sortIdx] = sort(R_locs_tmp,'descend');
R_angle_tmp = R_angle_tmp(sortIdx);
R_locs = mean(R_locs_tmp(1:round(0.1*length(R_locs_tmp))));

angle_range = [-180, 180];
angularBinNum = 36;
angularBins = linspace(angle_range(1),angle_range(2),angularBinNum)';
angularBinSize = angularBins(2)-angularBins(1);
angularBins = [angularBins];
angularRadiusMean = zeros(size(angularBins));

for ii = 1:angularBinNum
    if ii< angularBinNum-1
        angularDisplMean_ind = find(R_angle_tmp >= angularBins(ii) & R_angle_tmp < angularBins(ii+1));
    else
        angularDisplMean_ind = find((R_angle_tmp >= angularBins(1) & R_angle_tmp < angularBins(2)));
    end
    angularRadius_curr = R_locs_tmp(angularDisplMean_ind);
    angularRadiusMean_curr = mean(angularRadius_curr(1:round(0.1*length(angularRadius_curr))));
    
    angularRadiusMean(ii) = angularRadiusMean_curr;
end
Fig = figure('Position',[460,290,1300,375]); subplot(131), set(gca,'Color','k');
polarplot(angularBins*pi/180,angularRadiusMean,'LineWidth',2); 
thetatickformat('degrees')
title('Radius of the locs from the out layer'); 
subtitle('At different angle');


subplot(132), axis image; hold on;
th = 0:pi/50:2*pi;
xunit = R_locs * cos(th) + X_locs;
yunit = R_locs * sin(th) + Y_locs;
plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
hold on
scatter(SM_est_filtered(:,2),SM_est_filtered(:,3),2,'w','.');
title('Fitted circle and locs')
set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
xlim([X_locs-1.2*R_locs,X_locs+1.2*R_locs]); ylim([Y_locs-1.2*R_locs,Y_locs+1.2*R_locs]);
hold on,
plot([min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600),min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600)+500], [Y_locs-1.2*R_locs+100, Y_locs-1.2*R_locs+100],'w','LineWidth', 4);
% text(min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600)+250,Y_locs-1.2*R_locs+120,'500 nm','Color','w','FontSize',10,'FontWeight','bold','Interpreter','latex', 'HorizontalAlignment', 'center');
scatter(X_locs,Y_locs,20,'r','o','filled');


subplot(133),
SM_angles = atan2(SM_est_filtered(:,3)-X_locs,SM_est_filtered(:,2)-Y_locs);
polarhistogram(SM_angles, 54); title('Localization Density');
thetatickformat('degrees')


imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_initial_fit_circle_all_locs');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

clear th xunit yunit i ii angularDisplMean_ind angularRadiusMean_curr
clear angularRadiusMean angle_range angularBinNum angularBins angularBinSize
clear R_locs_tmp R_angle_tmp sortIdx angularRadius_curr SM_angles

%% Classify the loc into molecules
loc_list = SM_est_filtered(:,1:4);
loc_list(:,2:3) = loc_list(:,2:3);
frameN = max(loc_list(:,1));
threshold = 200;
check_frames = 2;

[loc_list, mol_ind] = classifyDiffuse(loc_list,threshold,check_frames,frameN);
%% Sort the molecules 
% mol_list. 1 - mol number; 2,3 - mean x and y; 4- mean brightness
% 5 - first frame; 6 - last frame;  7 - total lasting frames; 
% 8,9 - first x/y
molecule_list = zeros(max(loc_list(:,mol_ind)), 9);
molecule_list(:,1)=[1:max(loc_list(:,mol_ind))]'; 
for k = 1:max(loc_list(:,1))
    locs_list_tmp = loc_list(loc_list(:,5)==k,:);
    if ~isempty(locs_list_tmp) 
    firstFrameNum = min(locs_list_tmp(:,1),[],'all');
    lastFrameNum = max(locs_list_tmp(:,1),[],'all');
    trajFrameNum = lastFrameNum-firstFrameNum+1;
    meanX = mean(locs_list_tmp(:,2),'all');
    meanY = mean(locs_list_tmp(:,3),'all');
    meanS = mean(locs_list_tmp(:,4),'all');
    FirstX = locs_list_tmp(locs_list_tmp(:,1)==firstFrameNum,2);
    FirstY = locs_list_tmp(locs_list_tmp(:,1)==firstFrameNum,3);
    molecule_list(k,:)=[k, meanX, meanY, meanS, firstFrameNum, lastFrameNum, trajFrameNum, FirstX, FirstY];
    else 
        continue
    end
end
clear locs_list_tmp firstFrameNum lastFrameNum trajFrameNum meanX meanY meanS FirstX FirstY

%% brightness filter on molecules
br_threshold = 300;
molecule_list_br_filter_ind =  molecule_list(molecule_list(:,4)>br_threshold,1);
molecule_list = molecule_list(molecule_list_br_filter_ind,:);
clear molecule_list_br_filter_ind
%% Filter out the center
% R_coeff = 0.6;
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_coeff*R_locs)^2);
% molecule_list = molecule_list(molecule_list_R_filter_ind,:);
InterfaceThickness = 300; 
molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_locs-InterfaceThickness)^2);
molecule_list = molecule_list(molecule_list_R_filter_ind,:);
%% Second basic circle fit filter - for circle fit only!! Filter the centers
SM_circFit = molecule_list;
SM_circFit((SM_circFit(:,8)-X_locs).^2+(SM_circFit(:,9)-Y_locs).^2<(0.95*R_locs)^2,:) = [];
% SM_circFit(SM_circFit(:,3)<0,:) = [];
%% Secondary circular fit
% In fitting the curve, v2 updated using the nm scale coordinate (instead
% of pix) so the fitting circular is slightly different.
[R_mol,XC_mol,YC_mol] = circfit(SM_circFit(:,8),SM_circFit(:,9));
Fig = figure('Position',[600,200,800,375]);
subplot(121)
axis image; hold on;
th = 0:pi/50:2*pi;
xunit = R_mol * cos(th) + XC_mol;
yunit = R_mol * sin(th) + YC_mol;
plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
hold on
scatter(molecule_list(:,8),molecule_list(:,9),2,'w','.');
scatter(SM_circFit(:,8),SM_circFit(:,9),5,'r','.');
title('Fitted circle');subtitle('with trajactory start points');
set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
hold on,
plot([min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600),min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+500], [YC_mol-1.2*R_mol+100, YC_mol-1.2*R_mol+100],'w','LineWidth', 4);
scatter(XC_mol,YC_mol,20,'r','o','filled');

% text(min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+250,YC_mol-1.2*R_mol+220,'500 nm','Color','w','FontSize',10,'FontWeight','bold','Interpreter','latex','HorizontalAlignment', 'center');
% scatter(molecule_list(:,8),molecule_list(:,9),5,'r','.'); 

% Plot the molecular density across directions
molecule_list_angles = atan2(molecule_list(:,9)-YC_mol,molecule_list(:,8)-XC_mol);
subplot(122);
polarhistogram(molecule_list_angles, 54); title('Molecular Density');
thetatickformat('degrees')

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_secondary_fit_circle_all_locs');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

clear molecule_list_angles



%% Select dataset and interfaces (filters)
if dataNum == 2
    molecule_numel_cases = {molecule_list};
%     molecule_numel_cases(molecule_numel_cases(:,1)<4200,:) = [];
    XC_cases = XC_mol;
    YC_cases = YC_mol;
    molecule_numel_all = molecule_list;
    ImageNamingSuffix = {''};
    ArcColor = [199 160 151]/255; % red
elseif dataNum == 8
  
    switch fovNum
    case 1
        search_angle_arc_dil_dens = [-40,  140]; %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [170, 180; -180, -60];
        LineVec1_dens = [0, -1]; LineVec2_dens = [-1,0]; LineVec3_dens = [-1,-1]; % data 8 FoV 1
        LineVec1_dil = [1, 0]; LineVec2_dil = [0,1]; LineVec3_dil = [1,1]; % data 8 FoV 1
    case 2
        search_angle_arc_dil_dens =  [-20, 180; -180, -170]; %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [-150,  -30];
        LineVec1_dens = [1, -2]; LineVec2_dens = [0,-1]; LineVec3_dens = [-1,-2]; % data 8 FoV 2
        LineVec1_dil = [1, 2]; LineVec2_dil = [-1,2]; LineVec3_dil = [1,0];
    case 3
        search_angle_arc_dil_dens =  [-20,  165]; %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [-170, -35];
        LineVec1_dens = [0, -1]; LineVec2_dens = [1,-5]; LineVec3_dens = [-1,-1]; % data 8 FoV 3
        LineVec1_dil = [1, 0]; LineVec2_dil = [0,1]; LineVec3_dil = [1,1]; 
    case 4
        search_angle_arc_dil_dens = [150, 180; -180, -50]; %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [-25,  115];
        LineVec1_dens = [1, 0]; LineVec2_dens = [0,1]; LineVec3_dens = [1,1]; % data 8 FoV 4
        LineVec1_dil = [-1, -1]; LineVec2_dil = [0,-1]; LineVec3_dil = [-2,1]; 
    case 5
        search_angle_arc_dil_dens = [10,  160]; %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [-175, -30];
        LineVec1_dens = [0, -1]; LineVec2_dens = [1,-5]; LineVec3_dens = [-1,-1]; % data 8 FoV 5
        LineVec1_dil = [-1, 1]; LineVec2_dil = [0,1]; LineVec3_dil = [1,1]; 
    case 6
        search_angle_arc_dil_dens = [-60,  150];  %$ dont set a bound right at 180 when you have to separate there
        search_angle_arc_dens_dens = [170, 180; -180, -80];
        LineVec1_dens = [0, -1]; LineVec2_dens = [-1,0]; LineVec3_dens = [-1,-1]; % data 8 FoV 6
        LineVec1_dil = [1, 1]; LineVec2_dil = [0,1]; LineVec3_dil = [1,0];
    end
        

    Fig = figure('Position',[460,290,1300,375]); 
    subplot(131), hold on;
    axis image; 
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
    scatter(molecule_list(:,8),molecule_list(:,9),2,'w','.');
    scatter(SM_circFit(:,8),SM_circFit(:,9),5,'r','.');
    title('Segmenting the interfaces');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    plot([min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600),min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+500], [YC_mol-1.2*R_mol+100, YC_mol-1.2*R_mol+100],'w','LineWidth', 4);


    dilColor = [0.9290 0.6940 0.1250];
    denColor = [0.4940 0.1840 0.5560];
    for i = 1:size(search_angle_arc_dil_dens,1)
        for j = 1:size(search_angle_arc_dil_dens,2)
            if search_angle_arc_dil_dens(i,j)==180 || search_angle_arc_dil_dens(i,j)==-180
                continue
            else 
                line([XC_mol,XC_mol+(R_mol+300)*cosd(search_angle_arc_dil_dens(i,j))],...
                    [YC_mol,YC_mol+(R_mol+300)*sind(search_angle_arc_dil_dens(i,j))],...
                    'LineWidth',2,'color',dilColor)
            end
        end
    end
    for i = 1:size(search_angle_arc_dens_dens,1)
        for j = 1:size(search_angle_arc_dens_dens,2)
            if search_angle_arc_dens_dens(i,j)==180 || search_angle_arc_dens_dens(i,j)==-180
                continue
            else 
                line([XC_mol,XC_mol+(R_mol+300)*cosd(search_angle_arc_dens_dens(i,j))],...
                    [YC_mol,YC_mol+(R_mol+300)*sind(search_angle_arc_dens_dens(i,j))],...
                    'LineWidth',2,'color',denColor)
            end
        end
    end

    scatter(XC_mol,YC_mol,20,'r','o','filled');

    % %%% classify
    mol_polar_angle = atan2d(molecule_list(:,9)-YC_mol,molecule_list(:,8)-XC_mol);
    mol_dil_dens = [];
    mol_dens_dens = [];
    for i = 1:size(search_angle_arc_dil_dens,1)
        ind = find(mol_polar_angle>search_angle_arc_dil_dens(i,1) & mol_polar_angle<search_angle_arc_dil_dens(i,2));
        mol_dil_dens = [mol_dil_dens;molecule_list(ind,:)];
    end
    for i = 1:size(search_angle_arc_dens_dens,1)
        ind = find(mol_polar_angle>search_angle_arc_dens_dens(i,1) & mol_polar_angle<search_angle_arc_dens_dens(i,2));
        mol_dens_dens = [mol_dens_dens;molecule_list(ind,:)];
    end
    subplot(132), hold on;
    axis image; 
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    scatter(mol_dil_dens(:,8),mol_dil_dens(:,9),2,'w','.');
    title('dilute - dense interface'); 
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');
    plot([min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600),min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+500], [YC_mol-1.2*R_mol+100, YC_mol-1.2*R_mol+100],'w','LineWidth', 4);

    subplot(133), hold on;
    axis image; 
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    scatter(mol_dens_dens(:,8),mol_dens_dens(:,9),2,'w','.');
    title('dense - dense interface');  
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');
    plot([min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600),min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+500], [YC_mol-1.2*R_mol+100, YC_mol-1.2*R_mol+100],'w','LineWidth', 4);
    sgtitle('Using trajectory initial points')
    %%%

    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_segmenting_interfaces');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

    vec_normed = @(vec) vec./norm(vec);

    %%% dil-dens interface
    SM_circFit = mol_dil_dens;
    

    LineVec1 = vec_normed(LineVec1_dil); LineVec2 = vec_normed(LineVec2_dil); LineVec3 = vec_normed(LineVec3_dil);
    Proj1 = sort(SM_circFit(:,2:3)*LineVec1','descend');
    Proj2 = sort(SM_circFit(:,2:3)*LineVec2','descend'); 
    Proj3 = sort(SM_circFit(:,2:3)*LineVec3','descend');
    Point1 = mean(Proj1(1:round(0.01*length(Proj1))))*LineVec1;
    Point2 = mean(Proj2(1:round(0.01*length(Proj2))))*LineVec2;
    Point3 = mean(Proj3(1:round(0.01*length(Proj3))))*LineVec3;
    ThreePoints = [Point1;Point2;Point3];
    [R_tmp,XC_tmp,YC_tmp] = circfit(ThreePoints(:,1),ThreePoints(:,2));
    clear LineVec1 LineVec2 LineVec3 Proj1 Proj2 Proj3 Points1 Points2 Points3 ThreePoints
    
    SM_circFit((SM_circFit(:,8)-XC_tmp).^2+(SM_circFit(:,9)-YC_tmp).^2<(0.98*R_tmp)^2,:) = [];
    [R_dil,XC_dil,YC_dil] = circfit(SM_circFit(:,2),SM_circFit(:,3));
    clear R_tmp XC_tmp YC_tmp

    Fig = figure('Position',[600,200,800,375]);
    subplot(121)
    axis image; hold on;
    th = 0:pi/50:2*pi;
    xunit = R_dil * cos(th) + XC_dil;
    yunit = R_dil * sin(th) + YC_dil;
    plot(xunit, yunit,'Color',dilColor,'LineStyle','-','LineWidth',2);
    xlim([XC_dil-1.2*R_dil,XC_dil+1.2*R_dil]); ylim([YC_dil-1.2*R_dil,YC_dil+1.2*R_dil]);
    hold on
    scatter(mol_dil_dens(:,8),mol_dil_dens(:,9),2,'w','.');
    scatter(SM_circFit(:,8),SM_circFit(:,9),5,dilColor,'.');
    title('Fitted dilute-dense interface');subtitle('with trajactory start points');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');
    scatter(XC_dil,YC_dil,20,dilColor,'o','filled');
    plot([min(XC_dil+0.72*R_dil,XC_dil+1.2*R_dil-600),min(XC_dil+0.72*R_dil,XC_dil+1.2*R_dil-600)+500], [YC_dil-1.2*R_dil+100, YC_dil-1.2*R_dil+100],'w','LineWidth', 4);


     %%% dil-dens interface
    SM_circFit = mol_dens_dens;

    LineVec1 = vec_normed(LineVec1_dens); LineVec2 = vec_normed(LineVec2_dens); LineVec3 = vec_normed(LineVec3_dens);
    Proj1 = sort(SM_circFit(:,2:3)*LineVec1','descend');
    Proj2 = sort(SM_circFit(:,2:3)*LineVec2','descend'); 
    Proj3 = sort(SM_circFit(:,2:3)*LineVec3','descend');
    Point1 = mean(Proj1(1:round(0.01*length(Proj1))))*LineVec1;
    Point2 = mean(Proj2(1:round(0.01*length(Proj2))))*LineVec2;
    Point3 = mean(Proj3(1:round(0.01*length(Proj3))))*LineVec3;
    ThreePoints = [Point1;Point2;Point3];
    [R_tmp,XC_tmp,YC_tmp] = circfit(ThreePoints(:,1),ThreePoints(:,2));
    clear LineVec1 LineVec2 LineVec3 Proj1 Proj2 Proj3 Points1 Points2 Points3 ThreePoints
    
    SM_circFit((SM_circFit(:,8)-XC_tmp).^2+(SM_circFit(:,9)-YC_tmp).^2<(0.98*R_tmp)^2,:) = [];
    [R_dens,XC_dens,YC_dens] = circfit(SM_circFit(:,2),SM_circFit(:,3));
    clear R_tmp XC_tmp YC_tmp
    
    subplot(122)
    axis image; hold on;
    th = 0:pi/50:2*pi;
    xunit = R_dens * cos(th) + XC_dens;
    yunit = R_dens * sin(th) + YC_dens;
    plot(xunit, yunit,'Color',denColor,'LineStyle','-','LineWidth',2);
    xlim([XC_dens-1.2*R_dens,XC_dens+1.2*R_dens]); ylim([YC_dens-1.2*R_dens,YC_dens+1.2*R_dens]);
    hold on
    scatter(mol_dens_dens(:,8),mol_dens_dens(:,9),2,'w','.');
    scatter(SM_circFit(:,8),SM_circFit(:,9),5,denColor,'.');
    title('Fitted dense-dense interface');subtitle('with trajactory start points');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');
    scatter(XC_dens,YC_dens,20,denColor,'o','filled');
    plot([min(XC_dens+0.72*R_dens,XC_dens+1.2*R_dens-600),min(XC_dens+0.72*R_dens,XC_dens+1.2*R_dens-600)+500], [YC_dens-1.2*R_dens+100, YC_dens-1.2*R_dens+100],'w','LineWidth', 4);
    

    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_sepInterfaces_fittingCircle');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%%%%%   Assign values
    molecule_numel_cases = {mol_dil_dens,mol_dens_dens};
    XC_cases = [XC_dil;XC_dens];
    YC_cases = [YC_dil;YC_dens];
    molecule_numel_all = [mol_dil_dens;mol_dens_dens];
    ImageNamingSuffix = {'_Dil_Dense_Interface','_Dense_Dense_Interface'};
    ArcColor = [186 167 123; 175 165 184]/255; 
end

clear th x_unit y_unit
%%
loc_list_sorted_all =[];
for interfaceInd = 1:length(XC_cases)
    molecule_numel = molecule_numel_cases{interfaceInd};
    XC_use = XC_cases(interfaceInd);
    YC_use = YC_cases(interfaceInd);
    ImgSuffix = ImageNamingSuffix{interfaceInd};
    DisplmColor = ArcColor(interfaceInd,:); 
    
    molecule_list_selection = [];
    for ind = 1: size(molecule_numel,1)
        molecule_list_selection_tmp = molecule_list(molecule_list(:,1)==molecule_numel(ind),:);
        molecule_list_selection = [molecule_list_selection; molecule_list_selection_tmp];
    end
    clear molecule_list_selection_tmp
%% 
loc_list_filter_ind=[];
for i = 1:length(molecule_numel)
    loc_list_filter_ind_tmp = find(loc_list(:,mol_ind)==molecule_numel(i));
    loc_list_filter_ind=[loc_list_filter_ind;loc_list_filter_ind_tmp];
end

loc_list_filtered = loc_list(loc_list_filter_ind,:);
[loc_list_sorted, dist_ind, mol_angle_ind, seq_jump_dist_ind, ...
    seq_jump_angle_ind, seq_jump_par_ind,seq_jump_perp_ind, total_dist] = ...
    sortDisplacements(loc_list_filtered, mol_ind, XC_mol, YC_mol);
%% jump histograms
jump_list = loc_list_sorted(:,seq_jump_dist_ind);
jump_list = jump_list(jump_list~=0);
jump_par_list = loc_list_sorted(:,seq_jump_par_ind);
jump_par_list = jump_par_list(jump_par_list~=0);
jump_perp_list = loc_list_sorted(:,seq_jump_perp_ind);
jump_perp_list = jump_perp_list(jump_perp_list~=0);
jump_angle_local = atan2d(jump_perp_list,jump_par_list); 

%% JUMP historgram filtering (not used now)

% % mol_R = sqrt((loc_list_sorted(:,2)-XC).^2+(loc_list_sorted(:,3)-YC).^2);
% % 
% % figure;
% % axis equal;
% % hold on;
% % scatter(loc_list_sorted(:,2),loc_list_sorted(:,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% % title('All locs')
% % xlim([XC-max(mol_R)-200,XC+max(mol_R)+200]);ylim([YC-max(mol_R)-200,YC+max(mol_R)+200]);
% % 
% % figure;h = histogram(mol_R,30);title('Hist of locs'); xlabel('R to the center');
% % V = h.Values;
% % E = h.BinEdges;
% % % Use islocalmax
% % [maxH,L] = max(V);
% % left = E(L);
% % right = E(L+1);
% % center = (left + right)/2;
% % hold on; plot([center,center],[0,maxH+50],'r','LineWidth',2);
% % 
% % % left1 = E(L);
% % % right1 = E(L+1);
% % % center1 = (left1 + right1)/2;
% % % hold on; plot([center1,center1],[0,maxH+50],'r','LineWidth',2);
% % 
% % 
% % % left2 = E(L+2);
% % % right2 = E(L+3);
% % % center2 = (left2 + right2)/2;
% % % plot([center2,center2],[0,maxH+50],'r','LineWidth',2);

% % % center = 750;
% % % %%% % Choose the range of filtering
% % % % (alternative jump historgram with filtering)
% % % RfilterRange = [center, max(mol_R)+1]; %[min(mol_R)-1, center],[center, max(mol_R)+1], [center, center2]
% % % destinationFilter = find(mol_R>=RfilterRange(1)&mol_R<RfilterRange(2));
% % % jump_list = loc_list_sorted(destinationFilter,seq_jump_dist_ind);
% % % jump_list = jump_list(jump_list~=0);
% % % jump_par_list = loc_list_sorted(destinationFilter,seq_jump_par_ind);
% % % jump_par_list = jump_par_list(jump_par_list~=0);
% % % jump_perp_list = loc_list_sorted(destinationFilter,seq_jump_perp_ind);
% % % jump_perp_list = jump_perp_list(jump_perp_list~=0);
% % % jump_angle_local = atan2d(jump_perp_list,jump_par_list); 
% % % 
% % % figure;
% % % axis equal;
% % % hold on;
% % % scatter(loc_list_sorted(destinationFilter,2),loc_list_sorted(destinationFilter,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% % % title('Filtered locs')
% % % xlim([XC-max(mol_R)-200,XC+max(mol_R)+200]);ylim([YC-max(mol_R)-200,YC+max(mol_R)+200]);

%% Radius jump filter
% (alternative jump historgram with filtering)

% mol_R = sqrt((loc_list_sorted(:,2)-XC_mol).^2+(loc_list_sorted(:,3)-YC_mol).^2);
% OutRadius = sort(mol_R,'descend');
% OutRadius = mean(OutRadius(1:round(0.1*size(OutRadius,1))));
% % RfilterRange = [OutRadius-200, OutRadius+200];
% RfilterRange = [0.8*OutRadius, 1.2*OutRadius];
% % RfilterRange = [0.5*OutRadius, 0.8*OutRadius];
% destinationFilter = find(mol_R>=RfilterRange(1)&mol_R<RfilterRange(2));
% jump_list = loc_list_sorted(destinationFilter,seq_jump_dist_ind);
% jump_list = jump_list(jump_list~=0);
% jump_par_list = loc_list_sorted(destinationFilter,seq_jump_par_ind);
% jump_par_list = jump_par_list(jump_par_list~=0);
% jump_perp_list = loc_list_sorted(destinationFilter,seq_jump_perp_ind);
% jump_perp_list = jump_perp_list(jump_perp_list~=0);
% jump_angle_local = atan2d(jump_perp_list,jump_par_list); 
% % % plot the locs with this filter
% figure('Color',[0,0,0],'Position',[460,290,375,435]);  set(gcf,'color','none')
% axis image; hold on;  title('locs for diffusion analysis');
% scatter(loc_list_sorted(destinationFilter,2),loc_list_sorted(destinationFilter,3),2,'k','.'); 

%% Burst time and traj dist and displacements
burstTime = (molecule_list_selection(:,7)-1)*10;
Fig = figure('Position',[100,100,900,800]); subplot(311);
histogram(burstTime,'BinWidth',10.0001);xlabel('Burst Time (ms)'), ylabel('count'); title(['Burst time among tracked molecules, threshold = ', num2str(threshold), 'nm']), xlim([0,250])
subplot(312); histogram(total_dist,'BinWidth',20); xlabel('Distances (nm)'), ylabel('count'); title(['Total travel distance of molecules, threshold = ', num2str(threshold), 'nm']), xlim([0,1500])
subplot(313); histogram(jump_list,'Normalization','pdf','BinWidth',5);  xlabel('Displacement per 10ms'), ylabel('pdf'); title(['Displacement per 10 ms of molecules, threshold = ', num2str(threshold), 'nm']), xlim([0,300])

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_BurstTime_TrajDistance_Displacements');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%% molecule number v.s. time
slideWin = 100;
samplePoints = 150;
sampleTime = linspace(0, frameN-slideWin, samplePoints);
SampleMolNum = zeros(size(sampleTime));
for index = 1 :length(sampleTime)
    ind = find(molecule_list_selection(:,5)>sampleTime(index) & molecule_list_selection(:,5)<=sampleTime(index)+slideWin);
    SampleMolNum(index) = numel(ind)/slideWin;
end

Fig = figure('Position',[300,100,1000,255]);
plot(sampleTime,SampleMolNum,'LineWidth',2,'Color',[0.3 0.3 0.3]); xlabel('Frames'), ylabel('# of molecules/100 frames');
xticks([0 5000 10000 15000 20000]); xticklabels([0 5000 10000 15000 20000]);

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_moleculeNum over time');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);
%% Color palette
% perpColor= [0.8, 0.3, 0.3]; parColor= [0.09,0.64,0.66];
% perpTextColor= [0.55, 0.25, 0.25]; parTextColor= [0.03,0.4,0.4];
% perpColor= [0.5, 0.2, 0.2]; parColor= [0.03,0.35,0.35];
% perpTextColor= [0.2, 0.08, 0.08]; parTextColor= [0.01,0.15,0.15];

% perpColor= [204, 106, 20]/255; parColor= [20, 118, 204]/255;
% perpColor= [230, 48, 72]/255; parColor= 0.98-perpColor;
% perpColor= [248, 124, 100]/255; parColor= [100, 224, 248]/255;
% perpTextColor= 1-(1-perpColor).^0.5; parTextColor= 1-(1-parColor).^0.5;
% perpTextColor= perpColor.^2; parTextColor= parColor.^2;

% perpColor= [235 123 64]/255; parColor= [65 132 234]/255;
% perpTextColor= [214 87 21]/255; parTextColor =[10 92 214]/255;


% perpColor= [216 91 101]/255; parColor= [56 174 173]/255;
% perpTextColor= [207 69 70]/255; parTextColor =[12 150 149]/255;
% DisplmColor = [43 42 51]/255;

perpColor= [207 69 70]/255; parColor= [12 150 149]/255;
perpTextColor= [147 10 28]/255; parTextColor =[10 83 84]/255;
% DisplmColor = [165 140 162]/255;
% DisplmColor = [149 165 204]/255;
% DisplmColor = [175 165 184]/255;
% DisplmColor = [160 140 150]/255;

% color_tmp = [linspace(parTextColor(1),perpTextColor(1),21)',...
%         linspace(parTextColor(2),perpTextColor(2),21)',...
%         linspace(parTextColor(3),perpTextColor(3),21)'];
% 
% DisplmColor = color_tmp(10,:)*0.8;
% figure; colorbar; colormap(color_tmp)
clear color_tmp
%% 1D gauss fit
binsNum = 40; binsLimit = [-200,200];
pd_perp = fitdist(jump_perp_list,'Normal');
pd_par = fitdist(jump_par_list,'Normal');
ci99_perp = paramci(pd_perp,'Alpha',.01);
ci99_par = paramci(pd_par,'Alpha',.01);
Fig = figure('Position',[100,100,900,350]); 
h1 = histfit_YQdefined(jump_perp_list,"normal",binsLimit,binsNum);  hold on;
h2 = histfit_YQdefined(jump_par_list,"normal",binsLimit,binsNum); 
set(h1(1),'facecolor',perpColor,'faceAlpha',0.25,'DisplayName','perpendicular'); set(h1(2),'color', perpTextColor)
set(h2(1),'facecolor',parColor,'faceAlpha',0.25,'DisplayName','parallel'); set(h2(2),'color', parTextColor)
legend([h1(1),h2(1)]);

xlabel('Displacement per 10ms'), ylabel('counts');   
title('Displacement per frame of molecules'); xlim([-300,300]);

createtextbox_gaussfit_annotate(Fig,pd_perp,pd_par,ci99_perp,ci99_par,perpColor,parColor, perpTextColor, parTextColor);
% createtextbox_gaussfit_annotate_new(Fig,pd_perp,pd_par,ci99_perp,ci99_par,perpColor,parColor, perpTextColor, parTextColor);

% ax = gca;
% set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_1DGaussFit');
% imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_1DGaussFit_simpleVer');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%% New scatter plot - mean std 
range = [-220,220];
binsNum = 51; binsLimit = [-250,250];


Fig = figure('Position',[700,100,450,450]);
s1 = axes('Parent',Fig,'Position',[0.21,0.02,0.78,0.78]); hold on; %axis image;

R_plot = R_mol;
Y_plot = -R_plot;
x_circ_max = 600;
x_circ = linspace(-x_circ_max,x_circ_max,200);
y_circ = sqrt(R_plot^2-x_circ.^2)+Y_plot;

% fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], 0.8*[1 1 1],'EdgeColor','none');


% theta_list = atan2d(jump_perp_list,jump_par_list);
% theta_list(theta_list<=0&theta_list>=-180)=-theta_list(theta_list<=0&theta_list>=-180);
% theta_list(theta_list<=180&theta_list>=90)=180-theta_list(theta_list<=180&theta_list>=90);
% 
% cMap_angle = [linspace(parColor(1),perpColor(1),255)',...
%         linspace(parColor(2),perpColor(2),255)',...
%         linspace(parColor(3),perpColor(3),255)'];
% 
% theta_dic = linspace(0, 90, size(cMap_angle,1));
% for ind = 1: 2400
%     [~, color_ind] = min(abs(theta_list(ind)-theta_dic));
%     scatter(jump_par_list(ind),jump_perp_list(ind), 6,cMap_angle(color_ind,:),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);  
% end

scatter(jump_par_list,jump_perp_list, 3,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);  
% % plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
% line(range,[0,0]);
% line([0,0],range);
xlim(range); ylim(range);  axis off; box off; 


x_pos = -175; y_pos = -175; scal_pos=175;

quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^1.1,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^1.1,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'k')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'visible','off')
clear x_pos y_pos scal_pos

line([pd_par.mu pd_par.mu], [pd_perp.mu-2*pd_perp.sigma pd_perp.mu+2*pd_perp.sigma],'LineWidth',2,'Color', perpColor.^1.1,'LineStyle',':');
line([pd_par.mu pd_par.mu], [pd_perp.mu-pd_perp.sigma pd_perp.mu+pd_perp.sigma],'LineWidth',4,'Color', perpColor.^1.1,'Marker','_','MarkerSize',10);

line([pd_par.mu-2*pd_par.sigma pd_par.mu+2*pd_par.sigma],[pd_perp.mu pd_perp.mu], 'LineWidth',2,'Color', parColor.^1.1,'LineStyle',':');
line([pd_par.mu-pd_par.sigma pd_par.mu+pd_par.sigma],[pd_perp.mu pd_perp.mu],'LineWidth',4,'Color', parColor.^1.1,'Marker','|','MarkerSize',10);

% scatter(0,0,20,"black","filled");

s_perp = axes('Parent',Fig,'Position',[0.06,0.02,0.15,0.78]);
% s_perp = subplot(5,5,[6 11 16 21]);
h1 = histfit_YQdefined(jump_perp_list,"normal",binsLimit,binsNum);  view([-90,90])
xlim(range); 
set(s_perp,'XAxisLocation','bottom','YAxisLocation','origin','TickDir','out')
set(s_perp.YAxis,'Visible','off')

s_par = axes('Parent',Fig,'Position',[0.21,0.8,0.78,0.15]);
% s_par = subplot(5,5,[2 3 4 5]);
h2 = histfit_YQdefined(jump_par_list,"normal",binsLimit,binsNum);  
xlim(range); 
set(s_par,'XAxisLocation','top','YAxisLocation','origin','TickDir','out')
set(s_par.YAxis,'Visible','off')

set(h1(1),'facecolor',perpColor,'faceAlpha',0.7,'DisplayName','perpendicular'); set(h1(2),'color', perpTextColor)
set(h2(1),'facecolor',parColor,'faceAlpha',0.7,'DisplayName','parallel'); set(h2(2),'color', parTextColor)

set(s1,'DataAspectRatio',[1 1 1]);
linkaxes([s1,s_perp,s_par],'x');
set(s1,'YLim',s1.XLim);

createtextbox_gaussfit_annotate_SepData(Fig,pd_perp,pd_par,ci99_perp,ci99_par,perpColor, parColor, perpTextColor, parTextColor)


imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_1DGaussFit_scatterProj_meanStd');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%%  New scatter plot - quartile 

Fig = figure('Position',[700,100,450,450]);
s1 = axes('Parent',Fig,'Position',[0.21,0.05,0.75,0.75]);
R_plot = R_mol;
Y_plot = -R_plot;
x_circ_max = 600;
x_circ = linspace(-x_circ_max,x_circ_max,200);
y_circ = sqrt(R_plot^2-x_circ.^2)+Y_plot;

hold on; %axis image;
% fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.85 0.85 0.85],'EdgeColor','none');

scatter(jump_par_list,jump_perp_list, 3,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);  box off; 
% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
% line(range,[0,0]);
% line([0,0],range);
xlim(range); ylim(range);  axis off;

x_pos = -175; y_pos = -175; scal_pos=175;

quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^1.1,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^1.1,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'k')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'visible','off')
% clear x_pos y_pos scal_pos


line([median(jump_par_list) median(jump_par_list)], [prctile(jump_perp_list,1) prctile(jump_perp_list,99)],'LineWidth',2,'Color', perpColor.^1.1,'LineStyle',':');
line([median(jump_par_list) median(jump_par_list)], [prctile(jump_perp_list,10) prctile(jump_perp_list,90)],'LineWidth',4,'Color', perpColor.^1.1,'Marker','_','MarkerSize',10);

line([prctile(jump_par_list,1) prctile(jump_par_list,99)],[median(jump_perp_list) median(jump_perp_list)], 'LineWidth',2,'Color', parColor.^1.1,'LineStyle',':');
line([prctile(jump_par_list,10) prctile(jump_par_list,90)],[median(jump_perp_list) median(jump_perp_list)],'LineWidth',4,'Color', parColor.^1.1,'Marker','|','MarkerSize',10);

% scatter(0,0,50,"black","filled",'LineWidth',4);

% scatter(0,0,20,"black","filled");

s_perp = axes('Parent',Fig,'Position',[0.06,0.05,0.15,0.75]);
h1 = histfit_YQdefined(jump_perp_list,"normal",binsLimit,binsNum);  view([-90,90])
xlim(range); 
set(gca,'XAxisLocation','bottom','YAxisLocation','origin','TickDir','out')
set(gca().YAxis,'Visible','off')

s_par = axes('Parent',Fig,'Position',[0.21,0.8,0.75,0.15]);
h2 = histfit_YQdefined(jump_par_list,"normal",binsLimit,binsNum);  
xlim(range);  
set(gca,'XAxisLocation','top','YAxisLocation','origin','TickDir','out')
set(gca().YAxis,'Visible','off')


set(h1(1),'facecolor',perpColor,'faceAlpha',0.7,'DisplayName','perpendicular'); set(h1(2),'color', perpTextColor)
set(h2(1),'facecolor',parColor,'faceAlpha',0.7,'DisplayName','parallel'); set(h2(2),'color', parTextColor)

set(s1,'DataAspectRatio',[1 1 1]);
linkaxes([s1,s_perp,s_par],'x');
set(s1,'YLim',s1.XLim);

createtextbox_gaussfit_annotate_SepData(Fig,pd_perp,pd_par,ci99_perp,ci99_par,perpColor, parColor, perpTextColor, parTextColor)

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_1DGaussFit_scatterProj_prctiles');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

%% Pol plot of the displacements and displacement angle
Fig = figure('Position',[600,100,323,375]); %sgtitle('Diffusion in local coordinates','Color','w')
% subplot(121); 
polarhistogram((90-jump_angle_local)/180*pi, 36); title('Displacement \theta histogram');
thetatickformat('degrees');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
% set(gca,'ThetaColor','w','RColor','k'); 
% ax = gca;
% set(ax,'Color','k','RColor','w','ThetaColor','w'); set(gcf,'Color','none');
% set(ax.Title,'Color','w');
% subplot(122);  
% polarscatter((jump_angle_local)/180*pi,jump_list,4,'o','filled','MarkerFaceAlpha',0.2,'MarkerEdgeColor',"none"); 
% % polarscatter((jump_angle_local)/180*pi,jump_list,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
% title('Displacement');

% exportgraphics(gcf,strcat(ResultsSavePath,'\data_', num2str(dataNum),'FoV_', num2str(fovNum),'polarplots.pdf'),'Resolution',600);
imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_Displm_Angles');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);



%% Displacement plot (all displacements center at 0) 
displImgSz = ceil(max(max(abs(jump_par_list),[],'all'),max(abs(jump_perp_list),[],'all'))/30)*30;
diplBinNum = 24;
jumpNum = size(jump_par_list,1);

cMap = [linspace(1,DisplmColor(1)^2,256)', linspace(1,DisplmColor(2)^2,256)', linspace(1,DisplmColor(3)^2,256)'];

% Step 1, make 2D bins
% range = [-displImgSz,displImgSz];
x_ind = linspace(range(1),range(2),diplBinNum);
y_ind = linspace(range(1),range(2),diplBinNum);
binSize = x_ind(2)-x_ind(1);
[X,Y] = meshgrid(x_ind,y_ind);
z_count = zeros(size(X));
for ii = 1:diplBinNum
    for jj = 1:diplBinNum
        count_ind = find(jump_par_list >= x_ind(ii)-binSize/2 & jump_par_list < x_ind(ii)+binSize/2 & ...
            jump_perp_list >= y_ind(jj)-binSize/2 & jump_perp_list < y_ind(jj)+binSize/2);
        count_curr = numel(count_ind);
        z_count(jj,ii) = count_curr;

    end
end
% z_pdf = z_count/(dataSz*binSize);


Fig = figure('Position',[111,200,1700,375]); subplot(131)
hold on,

% fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.85 0.85 0.85],'EdgeColor','none');

scatter(jump_par_list(1:min(4000,length(jump_list))),jump_perp_list(1:min(4000,length(jump_list))),3,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0); axis image; box off; 
xlim(range); ylim(range); 
% title('Scatter of the displacements'); 
% xlabel('${\Delta d _\parallel}$','interpreter','latex', 'FontSize',18); 
% ylabel('${\Delta d _\perp}$','interpreter','latex', 'FontSize',18);
% ax = gca;  
% set(gcf,'Color','none');
% set(ax,'Color','k','XColor','k','YColor','k'); 
% set(ax,'YTickLabel',[],'YTickLabel',[]);
% set(ax.Title,'Color','k');
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% ax.XAxis.LineWidth = 2;
% ax.YAxis.LineWidth = 2;
% ax.XAxis.Color = parColor*0.8;
% ax.YAxis.Color = perpColor*0.8;
% ax.XAxis.FontSize = 13;
% ax.YAxis.FontSize = 13;

% x_pos = -180; y_pos = -180; scal_pos=175;
quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^1.1,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^1.1,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'k')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'visible','off')
% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
scatter(0,0,50,"black","filled",'LineWidth',4);
% clear x_pos y_pos scal_pos

% ax1 = subplot(122);  
% ax2 = copyobj(ax1,Fig);
% hold on;
% imagesc(ax1,x_ind,y_ind,z_count), 
% fitted_gauss = @(x,y) gaussmf(x,[pd_par.std,pd_par.mu]).*gaussmf(y,[pd_perp.std,pd_perp.mu]);
% [C,h]=contour(ax2,x_ind,y_ind,fitted_gauss(X,Y),4,'LineWidth',1.5,"ShowText",'off'); 
% % [C,h]=contour(ax2,x_ind,y_ind,z_count,4,'LineWidth',1,"ShowText",'on'); 
% % clabel(C,h,'FontSize',10,'Color','r');
% hcb = colorbar;
% % % Set colormaps
% NewGray = hot(20);
% NewGray = NewGray(1:10,:);
% colormap(ax1,cMap)
% colormap(ax2,NewGray), 
% % Set all other properties of ax1 before moving on
% % Finally, link the axis properties and turn off axis #2.
% ax2.UserData = linkprop([ax1,ax2],...
%     {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
%     'ydir','xdir','xlim','ylim'}); % add more props as needed
% ax2.Visible = 'off';
% axis image, box off;
% set(gca,'YDir','normal');
% title("2D displacement historgram"),xlim(range); ylim(range); 
% xlabel('Parallel displacement (unit: nm)'); ylabel('Perpendicular displacement (unit: nm)')
% ax = gca;
% set(ax,'XColor','k','YColor','k');
% set(ax.Title,'Color','k');
% hcb.Title.String  = 'counts';
% hold on; 

subplot(132);  hold on
imagesc(x_ind,y_ind,z_count), 
colormap(cMap), hcb1 = colorbar;
axis image, box off;
set(gca,'YDir','normal');
title("2D displacement historgram"),xlim(range); ylim(range); 
hcb1.Title.String  = 'counts';

% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
scatter(0,0,50,"black","filled",'LineWidth',4);

% x_pos = -180; y_pos = -180; scal_pos=175;
scatter(0,0,50,"black","filled",'LineWidth',4);
quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^1.1,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^1.1,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'k')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'visible','off')
% clear x_pos y_pos scal_pos



ax1 = subplot(133);  
ax2 = copyobj(ax1,Fig);
hold on;
imagesc(ax1,x_ind,y_ind,z_count), 

% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
scatter(0,0,50,"black","filled",'LineWidth',4);

hold on,
% x_pos = -180; y_pos = -180; scal_pos=175;
scatter(0,0,50,"black","filled",'LineWidth',4);
quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^1.1,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^1.1,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'k')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^1.1,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'visible','off')
% clear x_pos y_pos scal_pos

NewGray = copper(7);
NewGray = NewGray(end-4:end-1,:);
NewGray = flipud(NewGray);
th = 0:pi/50:2*pi;
sigmaList = [0.5, 1, 1.5, 2];
for sigmaInd = 1:length(sigmaList)
    xunit = sigmaList(sigmaInd) * pd_par.sigma * cos(th) + pd_par.mu;
    yunit = sigmaList(sigmaInd) * pd_perp.sigma * sin(th) + pd_perp.mu;
    CircP(sigmaInd) = plot(xunit, yunit,'Color',NewGray(sigmaInd,:),'LineStyle','-','LineWidth',2,"DisplayName",strcat('\pm',num2str(sigmaList(sigmaInd)),"\sigma"));
end
legend(CircP(1:end),'Orientation','horizontal','NumColumns',2);


hcb1 = colorbar(ax1);
% % Set colormaps
colormap(ax1,cMap)
hcb1.Title.String  = 'counts';

% Set all other properties of ax1 before moving on
% Finally, link the axis properties and turn off axis #2.

hold on

ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
axis image, box off;
set(gca,'YDir','normal');
title("2D displacement historgram"),xlim(range); ylim(range); 
% xlabel('Parallel displacement (unit: nm)'); ylabel('Perpendicular displacement (unit: nm)')
% ax = gca;
% set(ax,'XColor','k','YColor','k');
% set(ax.Title,'Color','k');



% hcb2 = colorbar(ax2,'North','Position',[0.7833,0.8350,0.0825,0.0262]);
% 
% colormap(ax2,NewGray), 
% hcb2.Ticks  = (sigmaList-0.25)/2;
% hcb2.TickLabels  = [sigmaList(1)+"\sigma",sigmaList(2)+"\sigma",sigmaList(3)+"\sigma",sigmaList(4)+"\sigma"];
% hcb2.FontSize = 10;



imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_DisplmScatters4000_2DHist');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

MatFileName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_intermediate');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));

%% Displacement scatters dark

ImagN = min(4000,length(jump_list));
Fig = figure('Position',[111,200,400,375]); 
scatter(jump_par_list(1:ImagN),jump_perp_list(1:ImagN), 3,sqrt(DisplmColor),'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0); axis image; box off; 
xlim(range); ylim(range); 

hold on,
x_pos = -180; y_pos = -180; scal_pos=175;
scatter(0,0,50,"w","filled",'LineWidth',4);
quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^0.9,'MaxHeadSize',5)
quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^0.9,'MaxHeadSize',5)
line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'w')
text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^0.9,'HorizontalAlignment','left','VerticalAlignment','top');
text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^0.9,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
set(gca,'XTickLabel',[],'YTickLabel',[]);
box off
% set(ax,'Color','k','XColor','k','YColor','k'); 
% set(ax,'YTickLabel',[],'YTickLabel',[]);
% set(ax.Title,'Color','k');
% set(gca,'visible','off')
set(gca,'Color','k')
clear x_pos y_pos scal_pos


imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_DisplmScatters_DarkBkg_4000Samples');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

MatFileName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_intermediate');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));


%% Displacement scatters dark (colorcoded) (optional)

% ImagN = min(4000,length(jump_list));
% Fig = figure('Position',[111,200,400,375]);  axis image;
% 
% xlim(range); ylim(range); 
% 
% hold on,
% x_pos = -180; y_pos = -180; scal_pos=175;
% 
% quiver(x_pos-2,y_pos,50,0,0,'LineWidth',2,'Color', parColor.^0.9,'MaxHeadSize',5)
% quiver(x_pos,y_pos-2,0,50,0,'LineWidth',2,'Color', perpColor.^0.9,'MaxHeadSize',5)
% line([scal_pos-50, scal_pos],[-scal_pos,-scal_pos],'LineWidth',4,'Color', 'w')
% text(x_pos, y_pos-4,  'Parallel','FontWeight','bold','fontsize', 9,'Color', parColor.^0.9,'HorizontalAlignment','left','VerticalAlignment','top');
% text(x_pos-8, y_pos, 'Perpendicular','FontWeight','bold', 'fontsize', 9,'Color', perpColor.^0.9,'HorizontalAlignment','left','VerticalAlignment','bottom','Rotation',90);
% set(gca,'XTickLabel',[],'YTickLabel',[]);
% box off
% % set(ax,'Color','k','XColor','k','YColor','k'); 
% % set(ax,'YTickLabel',[],'YTickLabel',[]);
% % set(ax.Title,'Color','k');
% % set(gca,'visible','off')
% set(gca,'Color','k')
% clear x_pos y_pos scal_pos
% 
% 
% cMap_curr = [linspace(parColor(1),perpColor(1),255)',...
%         linspace(parColor(2),perpColor(2),255)',...
%         linspace(parColor(3),perpColor(3),255)'];
% theta_dic = linspace(0, 90, size(cMap_curr,1));
% JumpX = jump_par_list(1:ImagN);
% JumpY = jump_perp_list(1:ImagN);
% JumpAng = atan2d(JumpY,JumpX);
% JumpAng(JumpAng<=0&JumpAng>=-180)=-JumpAng(JumpAng<=0&JumpAng>=-180);
% JumpAng(JumpAng<=180&JumpAng>=90)=180-JumpAng(JumpAng<=180&JumpAng>=90);
% 
% 
% for i = 1:length(JumpX)
%     [~, color_ind] = min(abs(JumpAng(i)-theta_dic));
%     scatter(JumpX(i,:),JumpY(i,:), 3,'filled','MarkerFaceColor',cMap_curr(color_ind,:),'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0); 
% end
% 
% scatter(0,0,50,"w","filled",'LineWidth',4);
% 
% clear cMap_curr theta_dic JumpX JumpY JumpAng color_ind
% 
% imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_DisplmScatters_DarkBkg_4000Samples_colored');
% exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);
% 
% MatFileName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_intermediate');
% save(strcat(MatResultsSavePath,MatFileName,'.mat'));

%% integrate loc_list_sorted
loc_list_sorted_all = [loc_list_sorted_all;loc_list_sorted];

end

%% Filtered trajectories
% Traj plotting
% To identify the mol_ind and the traj length

plot_range = [XC_mol-1.2*R_mol, XC_mol+1.2*R_mol, YC_mol-1.2*R_mol, YC_mol+1.2*R_mol];

burstTimeFiltVal = 8;
burstTimeFiltValUp = 20;
traj_brightness_thresh = 600;

molecule_list_selection_traj = [];
% sortrows(molecule_numel_all,1);
for ind = 1: size(molecule_numel_all,1)
    molecule_list_selection_tmp = molecule_list(molecule_list(:,1)==molecule_numel_all(ind,1),:);
    molecule_list_selection_traj = [molecule_list_selection_traj; molecule_list_selection_tmp];
end

% Filter molecules by
traj_length_mol_filt = molecule_list_selection_traj;
traj_length_mol_filt(traj_length_mol_filt(:,4)<=traj_brightness_thresh,:)=[];
traj_length_mol_filt(traj_length_mol_filt(:,7)<=1,:)=[];
traj_length_mol_filt(traj_length_mol_filt(:,7)<burstTimeFiltVal+1,:)=[];

traj_length_mol_filt(traj_length_mol_filt(:,7)>burstTimeFiltValUp+1,:)=[];

mol_select = traj_length_mol_filt(:,1);
colorselect = [0 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330];
colorselect = [colorselect, 0.7*ones(size(colorselect,1),1)]; % adjust tranparency

prename = strcat("s> ",num2str(traj_brightness_thresh), " ",num2str(burstTimeFiltVal*10)," <T< ",num2str(burstTimeFiltValUp*10),"ms");


Fig = plotTraj(loc_list,mol_select,mol_ind,colorselect,prename,plot_range);

imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_FilteredTraj');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

Fig = plotTraj_angle(loc_list_sorted_all,mol_select,mol_ind,seq_jump_par_ind,seq_jump_perp_ind,perpColor,parColor,prename,plot_range);
imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_FilteredTraj_colorcoded');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

% plotTrajGrowth_colorCoded(loc_list,mol_select,mol_ind,dataNum,threshold,check_frames,ResultsSavePath,15,prename,plot_range);
MatFileName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_all');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));


close all
end
close all
%% Displacement plot (all displacements center at 0) (extended version)
% displImgSz = ceil(max(max(abs(jump_par_list),[],'all'),max(abs(jump_perp_list),[],'all'))/30)*30;
% diplBinNum = 15;
% jumpNum = size(jump_par_list,1);
% 
% cMap = 1-copper;
% 
% % Step 1, make 2D bins
% range = [-displImgSz,displImgSz];
% x_ind = linspace(range(1),range(2),diplBinNum);
% y_ind = linspace(range(1),range(2),diplBinNum);
% binSize = x_ind(2)-x_ind(1);
% [X,Y] = meshgrid(x_ind,y_ind);
% z_count = zeros(size(X));
% densityList = [jump_par_list,jump_perp_list, zeros(size(jump_perp_list,1),4)];
% for ii = 1:diplBinNum
%     for jj = 1:diplBinNum
%         count_ind = find(jump_par_list >= x_ind(ii)-binSize/2 & jump_par_list < x_ind(ii)+binSize/2 & ...
%             jump_perp_list >= y_ind(jj)-binSize/2 & jump_perp_list < y_ind(jj)+binSize/2);
%         count_curr = numel(count_ind);
%         z_count(jj,ii) = count_curr;
%         if ~isempty(count_ind)
%             densityList(count_ind,3) = count_curr;
%         else 
%             continue
%         end
%     end
% end
% % z_pdf = z_count/(dataSz*binSize);
% 
% angle_range = [-180, 180];
% angularBinNum = 24;
% angularBins = linspace(angle_range(1),angle_range(2),angularBinNum)';
% angularBinSize = angularBins(2)-angularBins(1);
% angularBins = [angularBins];
% 
% angularDisplMean = zeros(size(angularBins));
% 
% 
% 
% for ii = 1:angularBinNum
%     if ii< angularBinNum-1
%         angularDisplMean_ind = find(jump_angle_local >= angularBins(ii) & jump_angle_local < angularBins(ii+1));
%     else
%         angularDisplMean_ind = find((jump_angle_local >= angularBins(1) & jump_angle_local < angularBins(2)));
%     end
%     angularDisplMean_curr = mean(jump_list(angularDisplMean_ind));
%     
%     angularDisplMean(ii) = angularDisplMean_curr;
% end
% 
% 
% Fig = figure('Position',[111,200,1750,385]); subplot(141)
% scatter(jump_par_list,jump_perp_list, 1,cMap(end,:),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0); axis image; box off; 
% xlim(range); ylim(range); title('Scatter of the diffusion data'); xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
% ax = gca;  
% % set(gcf,'Color','none');
% % set(ax,'Color','k','XColor','k','YColor','k'); 
% % set(ax,'YTickLabel',[],'YTickLabel',[]);
% set(ax.Title,'Color','k');
% 
% 
% subplot(143)
% polarplot(angularBins*pi/180,angularDisplMean,'LineWidth',2); 
% title('Average jump distance along angles'); %thetalabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
% ax = gca;  
% % set(gcf,'Color','none');
% set(ax,'Color','w','ThetaColor','k','RColor','k'); 
% % set(ax,'ThetaTickLabel',[],'RTickLabel',[]);
% set(ax.Title,'Color','k');
% %%%
% ax1 = subplot(142);  
% ax2 = copyobj(ax1,Fig);
% hold on;
% imagesc(ax1,x_ind,y_ind,z_count), 
% [C,h]=contour(ax2,x_ind,y_ind,z_count,4,'LineWidth',1,"ShowText",'on'); 
% clabel(C,h,'FontSize',10,'Color','r'); colorbar
% % Set colormaps
% colormap(ax1,cMap)
% colormap(ax2,'hot'), 
% % Set all other properties of ax1 before moving on
% % Finally, link the axis properties and turn off axis #2.
% ax2.UserData = linkprop([ax1,ax2],...
%     {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
%     'ydir','xdir','xlim','ylim'}); % add more props as needed
% ax2.Visible = 'off';
% axis image, box off;
% set(gca,'YDir','normal');
% title("Contour of fitted 2D Gaussian"),xlim(range); ylim(range); 
% xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
% ax = gca;
% set(ax,'Color','k','XColor','k','YColor','k');
% set(ax.Title,'Color','k');
% hold on; MultiNum = 2;
% 
% 
% subplot(144); hold on,
% densityRange = linspace(min(densityList(:,3)), max(densityList(:,3)), size(cMap,1));
% for i = 1: size(densityList,1)
%     [~,color_indx] = min(abs(densityList(i,3)-densityRange));
%     scatter(densityList(i,1),densityList(i,2),1,cMap(color_indx,:),'filled')
% end
% axis image, box off;
% set(gca,'YDir','normal');
% title("Contour of fitted 2D Gaussian"),xlim(range); ylim(range); 
% xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')

%% Displacement images (not used in current version)
% this figure is good to look at and validate if the curve fitting is good
% in separating the local coordinates (parallel/perpendicular) nicely. The
% edge is dominated with parallel move displacements only because of the
% geometry effects at the edge (santa view effect). This is not a strong
% evidence for us to to layers of the distribution. But might be good to
% eliminate the most out layer diffusions.
% 
% displacements = [];
% % x, y, distance, angle, global angle, displacement brightness
% for i = 1: size(loc_list_sorted,1)
%     if loc_list_sorted(i,seq_jump_dist_ind)~=0
%         jump_curr = loc_list_sorted(i,seq_jump_dist_ind);
%         jump_par_curr = loc_list_sorted(i,seq_jump_par_ind);
%         jump_perp_curr = loc_list_sorted(i,seq_jump_perp_ind);
%         jump_angle_local_curr = atan2d(jump_perp_curr,jump_par_curr); 
%         jump_angle_curr = atan2d(loc_list_sorted(i,3)-loc_list_sorted(i-1,3),...
%             loc_list_sorted(i,2)-loc_list_sorted(i-1,2));
%         
%         x_curr = (loc_list_sorted(i,2)+loc_list_sorted(i-1,2))/2;
%         y_curr = (loc_list_sorted(i,3)+loc_list_sorted(i-1,3))/2;
%         d_br_curr = (loc_list_sorted(i,4)+loc_list_sorted(i-1,4))/2;
%     else 
%         continue
%     end
%     displacements = [displacements; x_curr, y_curr, jump_curr, jump_angle_local_curr, jump_angle_curr,d_br_curr];
% end
% 
% displacements = displacements(displacements(:,6)>400,:);
% 
% cmap = redbluecmap;
% cmap = cmap(2:10,:);
% newCmap = [cmap; flip(cmap,1);cmap; flip(cmap,1)];
% newCmap = imresize(newCmap, [255, 3]);  % original color map contain just 11 colors, this increase it to 64
% 
% phi_dic = linspace(-180,180,255);
% figure('Position',[100,100, 600, 600]); hold on;
% for i = 1:size(displacements,1)
%     [~,color_indx] = min(abs(displacements(i,4)-phi_dic));
%     quiver(displacements(i,1),displacements(i,2),...
%         (displacements(i,3)/threshold).*cosd(displacements(i,5)),(displacements(i,3)/threshold).*sind(displacements(i,5)),...
%         100, 'ShowArrowHead',"off",'color',newCmap(color_indx,:));
% end
% axis image; colormap(newCmap); caxis([-180,180]), colorbar;
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% ax = gca;
% ax.FontSize = 10; 
% axis off;
% set(gcf, 'InvertHardcopy', 'off');
% set(gcf,'Color',[0 0 0]);

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
