
   clear; close all;
%% v5 diffuse analysis does not load raw images
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 2;
fovNum = 1;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');
% ResultsSavePath =[fullpath,'\ResultsSave - 20240312'];
% if ~exist(ResultsSavePath, 'dir')
%    mkdir(ResultsSavePath)
% end
%%
FigResultsSavePath = strcat("C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures\ProjectReview\fig_save",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'\');
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

% 
% InterfaceThickness = 300; 
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_locs-InterfaceThickness)^2);
% molecule_list_interface = molecule_list(molecule_list_R_filter_ind,:);
% 
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2<=(R_locs-InterfaceThickness)^2);
% molecule_list_center = molecule_list(molecule_list_R_filter_ind,:);
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

clear molecule_list_angles



%% Select dataset and interfaces (filters)
if dataNum == 2

    
    InterfaceThickness = 300; 
    molecule_list_R_filter_ind =  find((molecule_list(:,8)-XC_mol).^2+(molecule_list(:,9)-XC_mol).^2>(R_mol-InterfaceThickness)^2);
    molecule_list_interface = molecule_list(molecule_list_R_filter_ind,:);
    
    molecule_list_R_filter_ind =  find((molecule_list(:,8)-XC_mol).^2+(molecule_list(:,9)-YC_mol).^2<=(R_mol-InterfaceThickness)^2);
    molecule_list_center = molecule_list(molecule_list_R_filter_ind,:);  


    AreaCenter = 2*pi*((R_mol-InterfaceThickness)/1000)^2;   
    AreaInterface = 2*pi*((R_mol)/1000)^2-AreaCenter; 
    
    ArcColor = [199 160 151]/255; % red
    CenterColor = [175 175 175]/255;
    ImageNamingSuffix = {''};

    Fig = figure('Position',[600,200,800,375]);
    subplot(121)
    axis image; hold on;
    th = 0:pi/50:2*pi;
    xunit = R_mol * cos(th) + XC_mol;
    yunit = R_mol * sin(th) + YC_mol;
    plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    hold on
    scatter(molecule_list_interface(:,8),molecule_list_interface(:,9),2,'w','.');
    title('Fitted circle');subtitle('with trajactory start points');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on,
    plot([min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600),min(XC_mol+0.72*R_mol,XC_mol+1.2*R_mol-600)+500], [YC_mol-1.2*R_mol+100, YC_mol-1.2*R_mol+100],'w','LineWidth', 4);
    scatter(XC_mol,YC_mol,20,'r','o','filled');

    subplot(122)
    axis image; hold on;
    th = 0:pi/50:2*pi;
    xunit = R_mol * cos(th) + XC_mol;
    yunit = R_mol * sin(th) + YC_mol;
    plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    hold on
    scatter(molecule_list_center(:,8),molecule_list_center(:,9),10,'w','.');
    title('Fitted circle');subtitle('with trajactory start points');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');

    
    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_secondary_fit_circle_all_locs');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

    %%% molecule number v.s. time
    slideWin = 100;
    samplePoints = 150;
    sampleTime = linspace(0, frameN-slideWin, samplePoints);
    SampleMolNum_interface = zeros(size(sampleTime));
    SampleMolNum_center = zeros(size(sampleTime));
    for index = 1 :length(sampleTime)
        ind_interface = find(molecule_list_interface(:,5)>sampleTime(index) & molecule_list_interface(:,5)<=sampleTime(index)+slideWin);
        SampleMolNum_interface(index) = numel(ind_interface)/(slideWin*AreaInterface);
        ind_center = find(molecule_list_center(:,5)>sampleTime(index) & molecule_list_center(:,5)<=sampleTime(index)+slideWin);
        SampleMolNum_center(index) = numel(ind_center)/(slideWin*AreaCenter);
    end
    
    Fig = figure('Position',[300,100,1000,255]); hold on;
    plot(sampleTime,SampleMolNum_interface,'LineWidth',2,'Color',ArcColor,'DisplayName',"Interface"); 
    plot(sampleTime,SampleMolNum_center,'LineWidth',2,'Color',CenterColor,'DisplayName',"Interior");
    xticks([0 5000 10000 15000 20000]); xticklabels([0 50 100 150 200]);
    ylim([0 inf]); xlabel("time (s)"); ylabel("# of molecules \cdot \mu m ^{-2} \cdot s^{-1}");
    legend();
    
    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImageNamingSuffix,'_moleculeNum over time');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);
     %%%%


    molecule_numel_cases = {molecule_list_interface; molecule_list_center};
    XC_cases = XC_mol;
    YC_cases = YC_mol;
    molecule_numel_all = molecule_list;

    

elseif dataNum == 8

    InterfaceThickness = 300; 
    molecule_list_R_filter_ind =  find((molecule_list(:,8)-XC_mol).^2+(molecule_list(:,9)-YC_mol).^2<=(R_mol-InterfaceThickness)^2);
    molecule_list_center = molecule_list(molecule_list_R_filter_ind,:);  


    molecule_list_R_filter_ind =  find((molecule_list(:,8)-XC_mol).^2+(molecule_list(:,9)-XC_mol).^2>(R_mol-InterfaceThickness)^2);
    molecule_list = molecule_list(molecule_list_R_filter_ind,:);
    

    AreaCenter = 2*pi*((R_mol-InterfaceThickness)/1000)^2;   
    AreaInterface = 2*pi*((R_mol)/1000)^2-AreaCenter; 
    
    CenterColor = [175 175 175]/255;
    ImageNamingSuffix = {''};

    Fig = figure('Position',[600,200,800,375]);
    subplot(121)
  
    subplot(122)
    axis image; hold on;
    th = 0:pi/50:2*pi;
    xunit = R_mol * cos(th) + XC_mol;
    yunit = R_mol * sin(th) + YC_mol;
    plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
    xlim([XC_mol-1.2*R_mol,XC_mol+1.2*R_mol]); ylim([YC_mol-1.2*R_mol,YC_mol+1.2*R_mol]);
    hold on
    scatter(molecule_list_center(:,8),molecule_list_center(:,9),2,'w','.');
    title('Fitted circle');subtitle('with trajactory start points');
    set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
    scatter(XC_mol,YC_mol,20,'r','o','filled');

    
    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),'_secondary_fit_circle_inter_interior');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);
  
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
    CenterColor = [175 175 175]/255;
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
    Angle_dil_dens = 0;
    Angle_dens_dens = 0;
    for i = 1:size(search_angle_arc_dil_dens,1)
        ind = find(mol_polar_angle>search_angle_arc_dil_dens(i,1) & mol_polar_angle<search_angle_arc_dil_dens(i,2));
        mol_dil_dens = [mol_dil_dens;molecule_list(ind,:)];
        Angle_dil_dens = Angle_dil_dens + abs(search_angle_arc_dil_dens(i,2)-search_angle_arc_dil_dens(i,1));
    end
    for i = 1:size(search_angle_arc_dens_dens,1)
        ind = find(mol_polar_angle>search_angle_arc_dens_dens(i,1) & mol_polar_angle<search_angle_arc_dens_dens(i,2));
        mol_dens_dens = [mol_dens_dens;molecule_list(ind,:)];
        Angle_dens_dens = Angle_dens_dens + abs(search_angle_arc_dens_dens(i,2)-search_angle_arc_dens_dens(i,1));
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




    %%% molecule number v.s. time
    slideWin = 100;
    samplePoints = 150;
    sampleTime = linspace(0, frameN-slideWin, samplePoints);
    SampleMolNum_dil = zeros(size(sampleTime));
    SampleMolNum_den = zeros(size(sampleTime));
    SampleMolNum_center = zeros(size(sampleTime));
    for index = 1 :length(sampleTime)
        ind_dil = find(mol_dil_dens(:,5)>sampleTime(index) & mol_dil_dens(:,5)<=sampleTime(index)+slideWin);
        SampleMolNum_dil(index) = numel(ind_dil)/(slideWin*AreaInterface*Angle_dil_dens/360);
        ind_den = find(mol_dens_dens(:,5)>sampleTime(index) & mol_dens_dens(:,5)<=sampleTime(index)+slideWin);
        SampleMolNum_den(index) = numel(ind_den)/(slideWin*AreaInterface*Angle_dens_dens/360);
        ind_center = find(molecule_list_center(:,5)>sampleTime(index) & molecule_list_center(:,5)<=sampleTime(index)+slideWin);
        SampleMolNum_center(index) = numel(ind_center)/(slideWin*AreaCenter);
    end
    
    Fig = figure('Position',[300,100,1000,255]); hold on;
    plot(sampleTime,SampleMolNum_dil,'LineWidth',2,'Color',dilColor,'DisplayName',"Dense-Dilute Interface"); 
    plot(sampleTime,SampleMolNum_den,'LineWidth',2,'Color',denColor,'DisplayName',"Dense-Dense Interface"); 
    plot(sampleTime,SampleMolNum_center,'LineWidth',2,'Color',CenterColor,'DisplayName',"Interior");
    xticks([0 5000 10000 15000 20000]); xticklabels([0 50 100 150 200]);
    ylim([0 inf]); xlabel("time (s)"); ylabel("# of molecules \cdot \mu m ^{-2} \cdot s^{-1}");
    legend();
    
    imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImageNamingSuffix,'_moleculeNum over time');
    exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);
     %%%%

%%%%%   Assign values
    molecule_numel_cases = {mol_dil_dens,mol_dens_dens, molecule_list_center};
    XC_cases = [XC_dil;XC_dens];
    YC_cases = [YC_dil;YC_dens];
    molecule_numel_all = [mol_dil_dens;mol_dens_dens; molecule_list_center];
    ImageNamingSuffix = {'_Dil_Dense_Interface','_Dense_Dense_Interface'};
    ArcColor = [186 167 123; 175 165 184]/255; 


    
% InterfaceThickness = 300; 
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_locs-InterfaceThickness)^2);
% molecule_list_interface = molecule_list(molecule_list_R_filter_ind,:);
% 
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2<=(R_locs-InterfaceThickness)^2);
% molecule_list_center = molecule_list(molecule_list_R_filter_ind,:);


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


prename = strcat("s over ",num2str(traj_brightness_thresh), " ",num2str(burstTimeFiltVal*10)," T ",num2str(burstTimeFiltValUp*10),"ms");

plotTrajGrowth_colorCoded(loc_list,mol_select,mol_ind,dataNum,threshold,check_frames,FigResultsSavePath,15,prename,plot_range);


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
