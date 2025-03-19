clear; close all; clc;
%%
load('C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\data\082423 YOPRO\48\48_est_FoV1.mat');
%%
%% Filtering
SM_est_save_all(:,1) = Angle_save(:,1);

loc_list = SM_est_save_all(:,1:4);
%loc_list(:,2:3) = loc_list(:,2:3)/58.5+45;
loc_list(loc_list(:,4)<500,:)=[];
%% localization accumation map
img_sz = 91; pix_sz = 58.5; bin_sz = 10;
grid = (-(img_sz-1)/2 : bin_sz/pix_sz : (img_sz-1)/2) *pix_sz;
[X,Y] = meshgrid(grid,grid);
I_loc = zeros(size(X));
for i = 1:length(grid)
    for j = 1:length(grid)
        locs_temp = [];
        locs_temp = find( loc_list(:,2)>=X(i,j)-bin_sz/2 & ...
                          loc_list(:,2)<X(i,j)+bin_sz/2 &...
                          loc_list(:,3)>=Y(i,j)-bin_sz/2 &...
                          loc_list(:,3)<Y(i,j)+bin_sz/2);
        I_loc(i,j) = numel(locs_temp);
    end
end


%%
x = loc_list(:,2);
y = loc_list(:,3);
xCenter = mean(loc_list(:,2));
yCenter = mean(loc_list(:,3));
radiusData = sqrt((x-xCenter).^2+(y-yCenter).^2);
radiusNorm = radiusData./25;
%%
% r_tmp = reshape(radiusNorm,1,[]);
% i_tmp = reshape(cropData,1,[]);
% [R_tmp, ind]=sort(r_tmp);
% I_tmp = i_tmp(ind);
% plot(R_tmp,I_tmp);
%%
Radius_sample = 0:0.01:1.5;
radius_ind =  (Radius_sample(2:end)+Radius_sample(1:end-1))/2;
LocNum = zeros(size(radius_ind));
for i = 1:length(Radius_sample)-1
    DataPc = find(radiusNorm>=Radius_sample(i) & radiusNorm<Radius_sample(i+1));
    LocNum(i) = numel(DataPc);
end
%%
% figure('Position',[100,300,300,450]); tiledlayout('flow'); hold on;
% scatter(loc_list(:,2),loc_list(:,3),4,'white','o','filled','MarkerFaceAlpha',0.3);
% axis image; 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% ax = gca;
% ax.FontSize = 10; 
% axis off;
% set(gcf, 'InvertHardcopy', 'off');
% %whitebg('k');
% set(gcf,'Color',[0 0 0]);
% line([55 55+17.0940],[15 15],'Color','w','LineWidth',4);
% text(59,19,'1um','Color','w','FontWeight','bold','FontSize',15);

figure('Position',[100,300,300,450]); tiledlayout('flow'); hold on;
nexttile; imagesc(I_loc); colormap gray; hcb = colorbar; axis image; box off; axis off; 
% line([60 60+17.0940],[75 75],'Color','w','LineWidth',4);
% text(61,79,'1um','Color','w','FontWeight','bold');

set(gcf,'Color',[0 0 0]);
set(get(hcb,'Title'),'String','locs','fontsize',10, 'Color','w'); 
hcb.TickLabelInterpreter='latex';
hcb.FontSize = 11;
hcb.Color = 'w';
nexttile; hold on; xlabel('Distance from the center'), ylabel('Localization numbers');
plot(radius_ind, LocNum,'DisplayName', 'Average','LineWidth',2,'Color',[0 0.4470 0.7410]);
set(gca,'FontSize',11,'Color','none','XColor','w', 'YColor','w');