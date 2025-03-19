clc; clear; close all;

d = 200; % nm, displacement (arc length)
sigma = 50; % nm, presumed estimated std of the diffusion
R = 1000; % nm, radius of the condensate
whitebg;
%% part 0 : circle on a 3D sphere;
centerInit = [0;0;R];
Ry = @(theta) [cosd(theta), 0, sind(theta);     
                0,  1,  0;
               -sind(theta), 0, cosd(theta)];
beta = linspace(-180,180,1000);
circInit = [d*sind(beta);d*cosd(beta);sqrt(R^2-d^2)*ones(size(beta))];

rotAngle = 130; 
centerRot = Ry(rotAngle)*centerInit;
circRot = Ry(rotAngle)*circInit;

[X,Y,Z] = sphere(25); %
figure('Position',[100,100,600,600]); 
hold on     
surface(R*X,R*Y,R*Z,'FaceColor',[0.8,0.8,0.8],'FaceAlpha',0.2,...
            'FaceLighting', 'gouraud','EdgeColor',[0.8,0.8,0.8]);
plot3(circRot(1,:),circRot(2,:),circRot(3,:),'y','LineWidth',3); 
scatter3(centerRot(1),centerRot(2),centerRot(3),150,'y.'); axis image; 
plot(circRot(1,:),circRot(2,:),'r--','LineWidth',1.5); 
scatter(centerRot(1),centerRot(2),150,'r.'); axis image; 
% set(gcf, 'color', 'none');   
% set(gca, 'color', 'none');
view(127.5,30); axis equal; box off; 
xlabel('x','Position',[1.4,0]*R); 
ylabel('y','Position',[0,1.4]*R); 
zlabel('z','Position',[0,0,1.4]*R); 
set(gca,'FontSize',10,'FontWeight','bold');
ax = gca;                        % gets the current axes
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.LineWidth = 1.5;
ax.XAxis.FirstCrossoverValue  = 0; % X crossover with Y axis
ax.XAxis.SecondCrossoverValue  = 0; % X crossover with Z axis
ax.YAxis.FirstCrossoverValue  = 0; % Y crossover with X axis
ax.YAxis.SecondCrossoverValue  = 0; % Y crossover with Z axis
ax.ZAxis.FirstCrossoverValue  = 0; % Z crossover with X axis
ax.ZAxis.SecondCrossoverValue = 0; % Z crossover with Y axis
%% Suppl part 0: save as gif
% rotAngle = linspace(0,179,100);
% 
% for ii=1:length(rotAngle)
%     centerRot = Ry(rotAngle(ii))*centerInit;
%     circRot = Ry(rotAngle(ii))*circInit;
% 
% 
%     figure('Position',[100,100,600,600]); 
%     hold on     
%     surface(R*X,R*Y,R*Z,'FaceColor',[0.8,0.8,0.8],'FaceAlpha',0.2,...
%                 'FaceLighting', 'gouraud','EdgeColor',[0.8,0.8,0.8]);
%     plot3(circRot(1,:),circRot(2,:),circRot(3,:),'y','LineWidth',3); 
%     scatter3(centerRot(1),centerRot(2),centerRot(3),150,'y.'); axis image; 
%     plot(circRot(1,:),circRot(2,:),'r--','LineWidth',1.5); 
%     scatter(centerRot(1),centerRot(2),150,'r.'); axis image; 
%     % set(gcf, 'color', 'none');   
%     % set(gca, 'color', 'none');
%     view(127.5,30); axis equal; box off; 
%     xlabel('x','Position',[1.4,0]*R); 
%     ylabel('y','Position',[0,1.4]*R); 
%     zlabel('z','Position',[0,0,1.4]*R); 
%     set(gca,'FontSize',10,'FontWeight','bold');
%     ax = gca;                        % gets the current axes
%     ax.FontSize = 12;
%     ax.FontWeight = 'bold';
%     ax.LineWidth = 1.5;
%     ax.XAxis.FirstCrossoverValue  = 0; % X crossover with Y axis
%     ax.XAxis.SecondCrossoverValue  = 0; % X crossover with Z axis
%     ax.YAxis.FirstCrossoverValue  = 0; % Y crossover with X axis
%     ax.YAxis.SecondCrossoverValue  = 0; % Y crossover with Z axis
%     ax.ZAxis.FirstCrossoverValue  = 0; % Z crossover with X axis
%     ax.ZAxis.SecondCrossoverValue = 0; % Z crossover with Y axis
% 
%     frame1 = getframe(gcf);
%     im1 = frame2im(frame1);
%     [A1,map1] = rgb2ind(im1,256);
%     if ii == 1    
%         imwrite(A1,map1,'SphereProjected - Circle.gif')
%     else
%         imwrite(A1,map1,'SphereProjected - Circle.gif','WriteMode','append','DelayTime',0.1)
%     end
%     close all;
% end
%% Part 1: 1D simulation of theta on diffusion
alpha = asind(d/R);

theta1 = alpha:0.01:(90-alpha);
theta2 = (90-alpha): 0.01: 90;
theta3 = 90: 0.01: (90+alpha);
theta4 = (90+alpha):0.01: (180-alpha);

func1 = @(theta) R*abs(sind(theta)-sind(theta-alpha));
func2 = @(theta) R*abs(sind(theta)-sind(theta+alpha));
func3 = @(theta) R*abs(1-sind(theta));

figure; hold on;
l1 = plot(theta1,func1(theta1),'y','LineWidth',2,'DisplayName','Inward');
plot(theta2,func1(theta2),'y','LineWidth',2);
plot(theta3,func2(theta3),'y','LineWidth',2);
plot(theta4,func2(theta4),'y','LineWidth',2);

l2 = plot(theta1,func2(theta1),'r','LineWidth',2,'DisplayName','Outward');
plot(theta2,func3(theta2),'r','LineWidth',2);
plot(theta3,func3(theta3),'r','LineWidth',2);
plot(theta4,func1(theta4),'r','LineWidth',2);

legend([l1,l2]);
title(['Displacement per 10 ms = ', num2str(d),' nm']);
xlabel('\theta'); ylabel('2D projected perpendicular displacement')

%% Part 2: 3D to 2D realization projection (circles)
% assume x-y-z coordinate while we only consider the isotropic diffusion on
% the surface of a condensate, and it will rotate along the y-axis towards
% x-positive, and form a \theta angle to z-positive.
%
% all diffusion surrounds a center, which is the origin of the projected
% coordinates. Here are some conversion after projection:
%
% y-positive/negative -> parallel positive/negative (this can be reversed,
% doesn't matter)
% x-positive -> perpendicular outward
% x-negative -> perpendicular inward
%
% For simplicity, I use the projected length to generate a circle with 
% close enough arc lengths 
centerInit = [0;0;R];
Ry = @(theta) [cosd(theta), 0, sind(theta);     
                0,  1,  0;
               -sind(theta), 0, cosd(theta)];
beta = linspace(-180,180,1000);
circInit = [d*sind(beta);d*cosd(beta);sqrt(R^2-d^2)*ones(size(beta))];

rotAngle = 75; 
centerRot = Ry(rotAngle)*centerInit;
circRot = Ry(rotAngle)*circInit;

figure; hold on; plot(R*sind(beta),R*cosd(beta),'w','LineStyle','--');
plot(circInit(2,:),circInit(1,:),'r','LineWidth',1); 
scatter(centerInit(2),centerInit(1),'r.'); axis image; 
xlim([-1000,1000]);ylim([-1000,1000]);
title('Projected sotropic diffusion projection at center');
subtitle(['with displacement of ', num2str(d),' nm']);
xlabel('Parallel diffusion'); ylabel('Perpendicular diffusion')

figure('Position',[800,100,433,910]);
tiledlayout(2,1,'TileSpacing','compact'); 
nexttile;
hold on; plot(R*sind(beta),R*cosd(beta),'w','LineStyle','--');
plot(circRot(2,:),circRot(1,:),'r','LineWidth',1); 
scatter(centerRot(2),centerRot(1),'r.'); axis image; 
set(gca,'YDir','reverse');
xlim([-1000,1000]);ylim([-1000,1000]); hold off;
title('Projected sotropic diffusion projection');
subtitle(['with displacement of ', num2str(d),' nm, \theta = ', ...
    num2str(rotAngle),'^\circ']);
xlabel('Parallel diffusion (nm) as a circle'); 
ylabel('Perpendicular diffusion (nm) as a circle')
nexttile
hold on; plot(circRot(2,:),circRot(1,:),'r','LineWidth',2); axis image; 
set(gca,'YDir','reverse');
scatter(centerRot(2),centerRot(1),'r','filled');
xlabel('Par'); ylabel('Perp')
%% Save Part 2 as gif
% rotAngle = linspace(0,89,100);
% for ii=1:length(rotAngle)
%     centerRot = Ry(rotAngle(ii))*centerInit;
%     circRot = Ry(rotAngle(ii))*circInit;
% 
% 
%     figure('Position',[800,100,433,910]);
%     tiledlayout(2,1,'TileSpacing','compact'); 
%     nexttile;
%     hold on; plot(R*sind(beta),R*cosd(beta),'w','LineStyle','--');
%     plot(circRot(2,:),circRot(1,:),'r','LineWidth',1); 
%     scatter(centerRot(2),centerRot(1),'r.'); axis image; 
%     set(gca,'YDir','reverse');
%     xlim([-1000,1000]);ylim([-1000,1000]); hold off;
%     title('Projected sotropic diffusion projection');
%     subtitle(['with displacement of ', num2str(d),' nm, \theta = ', ...
%         num2str(rotAngle(ii)),'^\circ']);
%     xlabel('Parallel diffusion (nm) as a circle'); 
%     ylabel('Perpendicular diffusion (nm) as a circle')
%     nexttile
%     hold on; plot(circRot(2,:),circRot(1,:),'r','LineWidth',2); axis image; 
%     set(gca,'YDir','reverse');
%     scatter(centerRot(2),centerRot(1),'r','filled');
%     xlabel('Par'); ylabel('Perp')
% 
%     frame1 = getframe(gcf);
%     im1 = frame2im(frame1);
%     [A1,map1] = rgb2ind(im1,256);
%     if ii == 1    
%         imwrite(A1,map1,'SphereIsotropicDiffusion - Circle.gif')
%     else
%         imwrite(A1,map1,'SphereIsotropicDiffusion - Circle.gif','WriteMode','append','DelayTime',0.1)
%     end
%     close all;
% end
%% Part 3: Use gaussian distributed displacement, then 2D gauss fit
% Here use a simple 2D gaussian on a projected plane to simulate
% For complete gaussian distribution on sphere, check: 
%       5-parameter Fisher–Bingham distribution
dataN = 1000;
xInit = normrnd(0,sigma,[1,1000]);
yInit = normrnd(0,sigma,[1,1000]);
zInit = sqrt(R^2-(xInit.^2+yInit.^2));

DiffInit = [xInit;yInit;zInit];
binsNum = 40; binsLimit = [-200,200];
%%
% the same codes in part 2
centerInit = [0;0;R];
Ry = @(theta) [cosd(theta), 0, sind(theta);     
                0,  1,  0;
               -sind(theta), 0, cosd(theta)];
beta = linspace(-180,180,1000);
rotAngle = 79; 
centerRot = Ry(rotAngle)*centerInit;
DiffRot = Ry(rotAngle)*DiffInit;
% New plotting
figure; hold on; plot(R*sind(beta),R*cosd(beta),'y','LineStyle','--');
scatter(DiffInit(2,:),DiffInit(1,:),6,'r.','LineWidth',1); 
scatter(centerInit(2),centerInit(1),300,'w.'); axis image; 
xlim([-1000,1000]);ylim([-1000,1000]);
title('Projected sotropic diffusion projection at center');
subtitle(['with displacement std of ', num2str(d),' nm']);
xlabel('Parallel diffusion'); ylabel('Perpendicular diffusion')


figure('Position',[800,100,433,910]);
tiledlayout(2,1,'TileSpacing','compact'); 
nexttile;
hold on; plot(R*sind(beta),R*cosd(beta),'y','LineStyle','--');
scatter(DiffRot(2,:),DiffRot(1,:),6,'r.'); 
scatter(centerRot(2),centerRot(1),300,'w.'); axis image; 
set(gca,'YDir','reverse');
xlim([-1000,1000]);ylim([-1000,1000]); hold off;
title('Projected sotropic diffusion projection');
subtitle(['with displacement of ', num2str(d),' nm, \theta = ', ...
    num2str(rotAngle),'^\circ']);
xlabel('Parallel diffusion circle'); ylabel('Perpendicular diffusion circle')

nexttile
hold on; scatter(DiffRot(2,:),DiffRot(1,:),6,'r.'); 
axis image; 
set(gca,'YDir','reverse');
scatter(centerRot(2),centerRot(1),300,'w.');
xlabel('Par'); ylabel('Perp')

% 2D gaussian fitting
jump_par_list = (DiffRot(2,:) - centerRot(2))';
jump_perp_list = (DiffRot(1,:) - centerRot(1))';

[mu_x, mu_y, sigma_x, sigma_y, angle, fiterr, rr, fig1, fig2] = diffusion2Dgauss(jump_par_list,jump_perp_list);  
%%% this gaussian fitting will need replaced soon.

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
%% suppl part 3: recreate part 1 with the fitted results
