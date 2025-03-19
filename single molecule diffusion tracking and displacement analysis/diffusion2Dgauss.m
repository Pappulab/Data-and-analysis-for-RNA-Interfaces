function [mu_x, mu_y, sigma_x, sigma_y, angle, fiterr, rr, fig1, fig2] = diffusion2Dgauss(x,y);




ImgSz = ceil(max(max(abs(x),[],'all'),max(abs(y),[],'all'))/30)*30;
binNum = 30;
dataSz = size(x,1);


% Step 1, make 2D bins
range = [-ImgSz,ImgSz];
x_ind = linspace(range(1),range(2),binNum);
y_ind = linspace(range(1),range(2),binNum+1);
binSize = x_ind(2)-x_ind(1);
[X,Y] = meshgrid(x_ind,y_ind);
z_count = zeros(size(X));
for ii = 1:binNum
    for jj = 1:binNum
        count_curr = numel(find(x >= x_ind(ii)-binSize/2 & x < x_ind(ii)+binSize/2 & ...
            y >= y_ind(jj)-binSize/2 & y < y_ind(jj)+binSize/2));
        z_count(jj,ii) = count_curr;
    end
end
z_pdf = z_count/(dataSz*binSize);
% Step 2, gaussian fit
 [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(X,Y,z_count);

 mu_x = fitresult(:,5);
 mu_y = fitresult(:,6);
 sigma_x = fitresult(:,3);
 sigma_y = fitresult(:,4);
 angle = fitresult(:,2);

  mu_x_err = fiterr(:,5);
 mu_y_err = fiterr(:,6);
 sigma_x_err = fiterr(:,3);
 sigma_y_err = fiterr(:,4);

 

%%
fig1 = figure('Position',[100,100,1200,300]); subplot(131), 
scatter(x,y,'w','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0); axis image; box off
xlim(range); ylim(range); title('Scatter of the diffusion data'); xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
ax = gca; set(gcf,'Color','none');
set(ax,'Color','k','XColor','w','YColor','w'); 
set(ax.Title,'Color','w');

subplot(132),
surf(x_ind,y_ind,z_count);  box off, view(-28,43); colormap hot; title('2D hist of diffusion')
set(gca,'YDir','normal'); xlim(range); ylim(range); zlim([0,max(z_count,[],'all')+1]); xlabel('Par'); ylabel('Perp');
xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,0.5,1],'Rotation',26)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1.6,9,1],'Rotation',-48)
ax = gca; set(gcf,'Color','none');
set(ax,'Color','k','XColor','w','YColor','w','ZColor','w'); 
set(ax.Title,'Color','w');

subplot(133), 
surf(x_ind,y_ind,zfit);  box off, view(-28,43);  colormap hot; title('Fit 2D Gaussian')
set(gca,'YDir','normal'); xlim(range); ylim(range); zlim([0,max(z_count,[],'all')+1]);xlabel('Par'); ylabel('Perp');
xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,0.5,1],'Rotation',26)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1.6,9,1],'Rotation',-48)
ax = gca; set(gcf,'Color','none');
set(ax,'Color','k','XColor','w','YColor','w','ZColor','w'); 
set(ax.Title,'Color','w');
% 
% figure('Position',[100,100,400,400]); plot(x_ind,sum(z_count,1));
% figure('Position',[100,100,400,400]); plot(y_ind,sum(z_count,2));
%%
fig2 = figure('Position',[100,100,900,400]); subplot(121), 
scatter(x,y,'w','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0); axis image; box off
xlim(range); ylim(range); title('Scatter of the diffusion data'); xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
ax = gca; set(gcf,'Color','none');
set(ax,'Color','k','XColor','w','YColor','w'); 
set(ax.Title,'Color','w');

% Create main axis and then copy it.  
ax1 = subplot(122);  
ax2 = copyobj(ax1,fig2);
% Plot on both axes; specify axis handles!
imagesc(ax1,x_ind,y_ind,z_count), 
[C,h]=contour(ax2,x_ind,y_ind,zfit,4,'LineWidth',1,"ShowText",'on'); 
clabel(C,h,'FontSize',10,'Color','white')
% Set colormaps
colormap(ax1,'gray')
colormap(ax2,'cool'), 
% Set all other properties of ax1 before moving on
% Finally, link the axis properties and turn off axis #2.
ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
axis image, box off;
set(gca,'YDir','normal');
title("Contour of fitted 2D Gaussian"),xlim(range); ylim(range); 
xlabel('Parallel diffusion (unit: nm)'); ylabel('Perpendicular diffusion (unit: nm)')
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w');
set(ax.Title,'Color','w');
hold on; MultiNum = 2;
% line([mu_x-plot_times*sigma_x,mu_x-plot_times*sigma_x],[range(1)/2, range(2)],'LineStyle','--','LineWidth',1.5, 'Color',0.6*[1,1,1])
% line([mu_x+plot_times*sigma_x,mu_x+plot_times*sigma_x],[range(1)/2, range(2)],'LineStyle','--','LineWidth',1.5, 'Color',0.6*[1,1,1])
% 
% line([range(1)/1.5, range(2)],[mu_y-plot_times*sigma_y,mu_y-plot_times*sigma_y],'LineStyle','--','LineWidth',1.5, 'Color',0.6*[1,1,1])
% line([range(1)/1.5, range(2)],[mu_y+plot_times*sigma_y,mu_y+plot_times*sigma_y],'LineStyle','--','LineWidth',1.5, 'Color',0.6*[1,1,1])

% For text alignment
sigmaRange_min = 30;
sigmaRange_max = 50;

text(mu_x,mu_y+(MultiNum+1)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),...
    [' $\sigma_{\parallel}$ = ',num2str(sigma_x,'%.3f')],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'Rotation',0, ...
    'FontSize',14,...
    'Color','w');

arrow1 = annotation('doublearrow',[0.67 0.81],...
    [0.78 0.78],'Color',[1 1 1],'LineWidth',2);
arrow1.Parent = gca;           % associate the arrow the the current axes
arrow1.X = [mu_x-MultiNum*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_x+MultiNum*min([max([sigma_x,sigmaRange_min]),sigmaRange_max])]; 
arrow1.Y = [mu_y+(MultiNum+0.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]), mu_y+(MultiNum+0.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max])];

text(mu_x+(MultiNum+1)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_y,...
    [' $\sigma_{\perp}$ = ',num2str(sigma_y,'%.3f')],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'Rotation',-90, ...
    'FontSize',14,...
    'Color','w');

arrow2 = annotation('doublearrow',[0.67 0.81],...
    [0.78 0.78],'Color',[1 1 1],'LineWidth',2);
arrow2.Parent = gca;           % associate the arrow the the current axes
arrow2.X = [mu_x+(MultiNum+0.5)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_x+(MultiNum+0.5)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max])]; 
arrow2.Y = [mu_y-MultiNum*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),mu_y+MultiNum*min([max([sigma_y,sigmaRange_min]),sigmaRange_max])];

text(mu_x,mu_y-(MultiNum+1.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),...
    ['[$\mu_{\parallel}$,$\mu_{\perp}$] = [',num2str(mu_x,'%.3f'),', ',num2str(mu_y,'%.3f'),']'],...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'Rotation',0, ...
    'FontSize',14,...
    'Color','w');
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w'); set(gcf,'Color','none');
set(ax.Title,'Color','w');
% 
% text(mu_x,mu_y+(MultiNum+1)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),...
%     [' $\sigma_{\parallel}$ = ',num2str(sigma_x,'%.2f'),'$\pm$',num2str(sigma_x_err,'%.2f')],...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'Rotation',0, ...
%     'FontSize',14);
% 
% arrow1 = annotation('doublearrow',[0.67 0.81],...
%     [0.78 0.78],'Color',[1 1 1],'LineWidth',2);
% arrow1.Parent = gca;           % associate the arrow the the current axes
% arrow1.X = [mu_x-MultiNum*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_x+MultiNum*min([max([sigma_x,sigmaRange_min]),sigmaRange_max])]; 
% arrow1.Y = [mu_y+(MultiNum+0.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]), mu_y+(MultiNum+0.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max])];
% 
% text(mu_x+(MultiNum+1)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_y,...
%     [' $\sigma_{\perp}$ = ',num2str(sigma_y,'%.2f'),'$\pm$',num2str(sigma_y_err,'%.2f')],...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'Rotation',-90, ...
%     'FontSize',14);
% 
% arrow2 = annotation('doublearrow',[0.67 0.81],...
%     [0.78 0.78],'Color',[1 1 1],'LineWidth',2);
% arrow2.Parent = gca;           % associate the arrow the the current axes
% arrow2.X = [mu_x+(MultiNum+0.5)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max]),mu_x+(MultiNum+0.5)*min([max([sigma_x,sigmaRange_min]),sigmaRange_max])]; 
% arrow2.Y = [mu_y-MultiNum*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),mu_y+MultiNum*min([max([sigma_y,sigmaRange_min]),sigmaRange_max])];
% 
% text(mu_x,mu_y-(MultiNum+1.5)*min([max([sigma_y,sigmaRange_min]),sigmaRange_max]),...
%     ['[$\mu_{\parallel}$,$\mu_{\perp}$] = [',num2str(mu_x,'%.2f'),'$\pm$',num2str(mu_x_err,'%.2f'),', ',num2str(mu_y,'%.2f'),'$\pm$',num2str(mu_y_err,'%.2f'),']'],...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'Rotation',0, ...
%     'FontSize',14);
end