clear; close all;
%% v5 diffuse analysis does not load raw images
%%
diffusionPath = 'C:\Users_NotBackedUp\YQ_data_local\Poly RNA condensates\processed_results\Diffusion Analysis';
dataNum = 2;
fovNum = 1;
% for i = 1:6
%     fovNum=i;
resultFolder = '\20240306 YoPro polyA polyAC _ stack 20000';
dataFolder = ''; %'\data\
fullpath = strcat(diffusionPath, resultFolder, dataFolder,'\data',num2str(dataNum),'_FoV',num2str(fovNum)');
% ResultsSavePath =[fullpath,'\ResultsSave - 20240312'];
% if ~exist(ResultsSavePath, 'dir')
%    mkdir(ResultsSavePath)
% end
%%
MatResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\mat_save_July2024",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'_center\');
if ~exist(MatResultsSavePath, 'dir')
   mkdir(MatResultsSavePath)
end

FigResultsSavePath = strcat("C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\fig_save_July2024",'\data',num2str(dataNum),'_FoV',num2str(fovNum),'_center\');
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
% R_coeff = 0.7;
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_coeff*R_locs)^2);
% molecule_list = molecule_list(molecule_list_R_filter_ind,:);
WallThickness = 300; 
molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2<(R_locs-WallThickness)^2);
molecule_list = molecule_list(molecule_list_R_filter_ind,:);
% InterfaceThickness = 300; 
% molecule_list_R_filter_ind =  find((molecule_list(:,8)-X_locs).^2+(molecule_list(:,9)-Y_locs).^2>(R_locs-InterfaceThickness)^2);
% molecule_list = molecule_list(molecule_list_R_filter_ind,:);


Fig = figure('Position',[460,290,323,375]);  set(gca,'Color','k');
axis image; hold on;
th = 0:pi/50:2*pi;
xunit = R_locs * cos(th) + X_locs;
yunit = R_locs * sin(th) + Y_locs;
plot(xunit, yunit,'Color',[195 0 23]/255,'LineStyle','-','LineWidth',2);
hold on
scatter(molecule_list(:,2),molecule_list(:,3),2,'w','.');
title('Centroid molecules')
set(gca,'Color','k'); set(gca,'XTickLabel',[],'YTickLabel',[]);
xlim([X_locs-1.2*R_locs,X_locs+1.2*R_locs]); ylim([Y_locs-1.2*R_locs,Y_locs+1.2*R_locs]);
hold on,
plot([min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600),min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600)+500], [Y_locs-1.2*R_locs+100, Y_locs-1.2*R_locs+100],'w','LineWidth', 4);
% text(min(X_locs+0.72*R_locs,X_locs+1.2*R_locs-600)+250,Y_locs-1.2*R_locs+120,'500 nm','Color','w','FontSize',10,'FontWeight','bold','Interpreter','latex', 'HorizontalAlignment', 'center');
scatter(X_locs,Y_locs,20,'r','o','filled');
%% Select dataset and interfaces (filters)
    molecule_numel_cases = {molecule_list};
%     molecule_numel_cases(molecule_numel_cases(:,1)<4200,:) = [];
    XC_cases = X_locs;
    YC_cases = Y_locs;
    molecule_numel_all = molecule_list;
    ImageNamingSuffix = {''};

clear th x_unit y_unit
%%
loc_list_sorted_all =[];
for interfaceInd = 1:length(XC_cases)
    molecule_numel = molecule_numel_cases{interfaceInd};
    XC_use = XC_cases(interfaceInd);
    YC_use = YC_cases(interfaceInd);
    ImgSuffix = ImageNamingSuffix{interfaceInd};
    
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
    sortDisplacements(loc_list_filtered, mol_ind, X_locs, Y_locs);
%% jump histograms
jump_list = loc_list_sorted(:,seq_jump_dist_ind);
jump_list = jump_list(jump_list~=0);
jump_par_list = loc_list_sorted(:,seq_jump_par_ind);
jump_par_list = jump_par_list(jump_par_list~=0);
jump_perp_list = loc_list_sorted(:,seq_jump_perp_ind);
jump_perp_list = jump_perp_list(jump_perp_list~=0);
jump_angle_local = atan2d(jump_perp_list,jump_par_list); 


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
plot(sampleTime,SampleMolNum,'LineWidth',2,'Color',[0.3 0.3 0.3]); xlabel('Frames'), ylabel('# of molecules');
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
DisplmColor = [175 175 175]/255;
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

R_plot = R_locs;
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

scatter(jump_par_list,jump_perp_list, 10,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);  
% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");
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
R_plot = R_locs;
Y_plot = -R_plot;
x_circ_max = 600;
x_circ = linspace(-x_circ_max,x_circ_max,200);
y_circ = sqrt(R_plot^2-x_circ.^2)+Y_plot;

hold on; %axis image;
% fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.85 0.85 0.85],'EdgeColor','none');

scatter(jump_par_list,jump_perp_list, 10,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0);  box off; 
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
polarhistogram((jump_angle_local)/180*pi, 36); title('Angles histogram');
thetatickformat('degrees');
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

scatter(jump_par_list(1:min(4000,length(jump_list))),jump_perp_list(1:min(4000,length(jump_list))),10,DisplmColor,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0); axis image; box off; 
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



imgName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_DisplmScatters4000_2DHist');
exportgraphics(Fig,strcat(FigResultsSavePath,imgName,'.pdf'),'ContentType','vector','Resolution',600);

MatFileName = strcat('data',num2str(dataNum),'_FoV',num2str(fovNum),ImgSuffix,'_intermediate');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));

%% Displacement scatters dark

ImagN = min(4000,length(jump_list));
Fig = figure('Position',[111,200,400,375]); 
scatter(jump_par_list(1:ImagN),jump_perp_list(1:ImagN), 10,sqrt(DisplmColor),'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0); axis image; box off; 
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


%% integrate loc_list_sorted
loc_list_sorted_all = [loc_list_sorted_all;loc_list_sorted];

end

%% Filtered trajectories
% Traj plotting
% To identify the mol_ind and the traj length

plot_range = [X_locs-1.2*R_locs, X_locs+1.2*R_locs, Y_locs-1.2*R_locs, Y_locs+1.2*R_locs];

burstTimeFiltVal = 8;
burstTimeFiltValUp = 20;
traj_brightness_thresh = 400;

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
% end
% close all

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
