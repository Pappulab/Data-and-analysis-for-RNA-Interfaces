clc;
clear;
% % dataXfovX = [pd_perp.mu, pd_perp.std, pd_par.mu, pd_par.std];
%% record data in here (archive on 03202024)
% % dataXfovX = [pd_perp_mu, pd_perp_std, pd_par_mu, pd_par_std];
% % poly-rA
% data2fov1 = [-9.88915, 60.0922, 0.472261, 71.2736]; % 21805 filtered locs -> 14780 displacement
% data2fov2 = [-7.96487, 60.3473, 0.346238, 69.7554]; % 10938 filtered locs -> 7498 displacement
% data2fov3 = [-9.45419, 61.1244, -0.064221, 70.6412]; % 10938 filtered locs -> 7498 displacement
% data2fov4 = [-8.80277, 60.0802 , -0.957763, 67.0644]; % 8833 filtered locs -> 5802 displacement
% data2fov5 = [-5.5243, 58.2734 , 0.2114, 67.2948]; % 23923 filtered locs -> 14249 displacement
% data2fov6 = [-6.5121, 58.6493 , -0.3906, 67.9123]; % 8833 filtered locs -> 5802 displacement
% 
% % poly-A/C - dense/dilute interface
% data8fov1_dilute = [-8.01757, 56.6587, 0.105483, 69.2333]; % 28134 filtered locs -> 17533 select locs  -> 11713 displacement
% data8fov2_dilute = [-8.39647, 55.2783, -0.215987, 68.8026]; % 42977 filtered locs -> 27476 select locs  -> 18613 displacement
% data8fov3_dilute = [-5.5103,  57.6529, 0.5742, 70.2485]; %  37917 filtered locs -> 21192 select locs  ->  13093 displacement
% data8fov4_dilute = [-5.7822,  56.4833, -0.7856, 67.4243]; % 19462 filtered locs -> 8835 select locs  -> 5156 displacement
% data8fov5_dilute = [-6.0164, 58.2513, -1.7599, 70.9407]; % 30205 filtered locs -> 13209 select locs  -> 8056 displacement
% data8fov6_dilute = [-9.8760, 55.8450, -0.2056, 69.1344]; % 31400 filtered locs -> 23461 select locs  -> 15560 displacement
% 
% % poly-A/C - dense/dense interface
% data8fov1_dense = [-8.61048, 56.8856, 2.40248, 65.6279]; % 28649 filtered locs -> 6690 select locs  -> 4337 displacement
% data8fov2_dense = [-8.14283, 55.5028, 0.948663, 68.0059]; % 21805 filtered locs -> 11132 select locs  -> 7552 displacement
% data8fov3_dense = [-6.2118, 56.5870, -0.3629, 67.1111]; % 37917 filtered locs -> 11345 select locs  -> 7143 displacement
% data8fov4_dense = [-5.5064, 57.9408, -1.0704, 67.7246]; % 19462 filtered locs -> 6551 select locs  -> 4057 displacement
% data8fov5_dense = [-6.1125, 58.0106, -0.6116, 69.1531]; % 30205 filtered locs -> 10652 select locs  -> 6669 displacement
% data8fov6_dense = [-9.7054, 54.2694, -0.5028, 64.3974]; % 31400 filtered locs -> 4262 select locs  -> 2766 displacement
% 

%% archives 20240401
% data_pA = [...
%     -10.4457185491371	59.9174518845295	0.491585183076329	71.4104894832901;... %data2 fov1
%     -8.66586499482549	60.0707976504147	0.379980394023473	69.7853514385214;... %data2 fov3
%     -10.1480887570246	60.7208268934696	0.138714505044242	70.8708720859766;... %data2 fov3
%     -9.08542830805685	59.8476228470827	-0.651203335940006	67.3755523814805;... %data2 fov4
%     -5.92457098862328	58.2801637936987	0.297309053716991	67.3352266755272;... %data2 fov5
%     -7.46805544320773	58.4858737383737	-0.387275736418306	68.2084344092763];   %data2 fov6
% data_pAC_dilute = [...
%     -8.58806929633013	56.6744669048224	0.223798258927945	69.7519527865479;... %data8 fov1 dil-dense
%     -8.59484173958744	55.1721386268426	-0.323063108528591	69.4353471723666;... %data8 fov2 dil-dense
%     -7.16834471083744	57.5603290284730	0.0598445975148858	71.4453030432343;... %data8 fov3 dil-dense
%     -7.30830739760044	56.3152576754586	-0.965319326324268	67.9770944314955;... %data8 fov4 dil-dense
%     -7.80623139941412	58.2080109081075	-1.72503615565139	72.2352209353512;... %data8 fov5 dil-dense
%     -9.96142877434737	55.8141503818743	-0.744919870173135	69.4163545203878];   %data8 fov6 dil-dense
% data_pAC_dense = [...
%     -8.39865510121827	57.4383385786471	0.932205765338908	66.5052299232206;... %data8 fov1 dense-dense
%     -8.36086784035328	55.2306953490672	0.993603258762359	69.2404833161449;... %data8 fov2 dense-dense
%     -7.58167812108896	55.9430755206063	0.100040663390059	68.8549303909478;... %data8 fov3 dense-dense
%     -7.02903012637729	57.7724999922259	0.405737379674338	68.6930629519161;... %data8 fov4 dense-dense
%     -7.34547109264620	57.8978189971687	0.152586626430636	70.3319599811618;... %data8 fov5 dense-dense
%     -9.32548329524272	54.3667656504501	0.311591435978006	65.5685022948159];   %data8 fov6 dense-dense

%%
data_pA = [...
    -11.8717536508830	59.2897434779593	0.847201953406734	71.8910639673408;... %data2 fov1
    -9.87719329527907	59.5968783466060	0.178698353220908	70.2147600699087;... %data2 fov2
    -11.5062626647334	60.0996643192910	0.172646598390231	71.2592075010529;... %data2 fov3
    -11.0461520727166	58.6345455900120	-0.546751168758521	67.9646338449999;... %data2 fov4
    -8.43296208954686	56.7206060642873	0.0582191940274565	67.9974244461784;... %data2 fov5
    -9.06426022328209	57.3926932670637	-0.785676082547558	68.7902279204900];   %data2 fov6
data_pAC_dilute = [...
    -9.64444268346021	55.6906581519866	0.140006568390980	70.1907874573259;... %data8 fov1 dil-dense
    -9.77491691707286	53.7019873801577	-0.373707382356111	69.7912119566707;... %data8 fov2 dil-dense
    -7.76587511481407	57.0689191876785	0.244736422561374	71.6405182895601;... %data8 fov3 dil-dense
    -8.33337385497666	55.8590194328317	-1.13029390227492	68.1295150545832;... %data8 fov4 dil-dense
    -8.42836104270591	57.6612501099518	-1.38586551836882	72.2498828511650;... %data8 fov5 dil-dense
    -11.1443983219478	54.7556796028291	-0.672628190820923	69.8865879664063];   %data8 fov6 dil-dense
data_pAC_dense = [...
    -9.95576258847996	56.0039174777834	0.672600092404191	66.8213097001604;... %data8 fov1 dense-dense
    -10.1985300896940	53.7659715499069	1.05397904463066	69.3970684625417;... %data8 fov2 dense-dense
    -8.45714288735698	55.3738974470255	0.142496893912633	68.9264943294344;... %data8 fov3 dense-dense
    -8.72929594686483	56.7679463678691	0.380633596148314	68.4752220030690;... %data8 fov4 dense-dense
    -8.21412544979082	57.1484446625007	-0.0918348370767480	70.0362552262432;... %data8 fov5 dense-dense
    -11.1229980322743	52.4209921246758	1.38805082212181	65.6756880928280];   %data8 fov6 dense-dense

%%
data_pA_center = [...
    6.69661571328870	64.2442218342126	-2.94579499623708	65.5584469605323;...
    6.43652365052650	63.6279293217139	1.68074518459540	66.0927843908356;...
    4.96478475065242	66.2865533442182	-1.95051332886047	66.0516293258239;...
    3.29063620738133	65.5866735227434	-2.59237208917513	62.6871308726359;...
    7.40378158474125	64.1907797080994	1.15998337747781	63.8052363072165;...
    4.54920003502465	62.9971068548169	1.06707715110934	64.0910497199180];   %data2 fov6
data_pAC_center = [...
    6.33390950817566	65.7804431257853	1.34615552238344	65.0924508306108;...
    3.74280110494122	65.4139435433709	0.642034289403216	65.4998886365753;...
    12.5681862240026	66.5521651183268	0.202867444102089	67.8684046382729;...
    10.8316435120228	63.2408539402345	-0.703158116349525	66.5362610897537;...
    12.9210202811805	67.7452656028704	-1.05959269385531	70.2211495605060;...
    2.86349736182797	64.2307282675513	-1.08854038602411	63.8574020089196];  

%%
% simData =[0; 35; 45; 55; 65; 75];
% simDataResults = [...
%      -0.1948   23.4026    0.4606   23.4110;... %static
%      2.4169   45.4284    0.2014   46.2370;...  %iso 35
%      3.7593   54.1552    0.7251   54.2738;...  %iso 45
%      4.2726   61.3313    0.1789   59.9501;...  %iso 55
%      5.5017   67.5373   -0.4018   68.7968;...  %iso 65
%      4.1236   73.2598   -0.5623   73.8766];    %iso 75

%% 4000 jumps per group (0320 data)
% simDataGroups =[0; 0; 35;35;35; 45;45;45; 55;55;55; 65;65;65; 75;75;75];
% simDataGroupsResults = [...
%      -0.4295   23.4496    0.4291   23.3129;... %static
%      0.4896   23.2410    0.0879   23.7718;... %static
%      ...
%      3.1106   45.5338    0.3752   46.3944;...  %iso 35
%      0.8316   45.4273    0.3900   45.6090;...  %iso 35
%      1.5452   45.9119   -0.7354   45.8865;...  %iso 35
%      ...
%      3.9577   54.3653    0.9841   54.0932;...  %iso 45
%      2.6931   53.2650   -0.2256   53.3906;...  %iso 45
%      2.9846   53.9467    0.3723   54.2595;...  %iso 45
%      ...
%      4.2456   61.6943   -0.1278   60.1065;...  %iso 55
%      4.1750   60.7749    0.0253   61.2117;...  %iso 55
%      2.6828   60.7500   -1.0491   62.0551;...  %iso 55
%      ...
%      5.7915   67.3532   -0.4916   68.2567;...  %iso 65
%      4.6038   68.8346    0.2200   69.6483;...  %iso 65
%      4.9698   68.1026   -0.0030   68.4435;...  %iso 65
%      ...
%      3.9979   72.6560    0.0242   73.8703;...  %iso 75
%      5.3560   75.4868   -1.3309   74.0779;...  %iso 75
%      3.5319   74.8857    3.4400   74.9433];    %iso 75

%% 4000 jumps per group (0328 data)
simDataGroups =[0;0; 0;0;0; 35;35;35; 45;45;45; 50;50;50;50;50; 55;55;55; ...65;65;65;
    65;65;65;65;65; 75;75;75];
simDataGroupsResults = [...
     -0.4295   23.4496    0.4291   23.3129;... %static
     0.4896   23.2410    0.0879   23.7718;... %static
     ...
     0.7009   22.8729    0.0505   23.3415;... %static
     0.1377   23.3845   -0.4889   23.7418;... %static
     0.0411   23.4858    0.2888   23.8098;... %static
     ...
     3.1106   45.5338    0.3752   46.3944;...  %iso 35
     0.8316   45.4273    0.3900   45.6090;...  %iso 35
     1.5452   45.9119   -0.7354   45.8865;...  %iso 35
     ...
     3.9577   54.3653    0.9841   54.0932;...  %iso 45
     2.6931   53.2650   -0.2256   53.3906;...  %iso 45
     2.9846   53.9467    0.3723   54.2595;...  %iso 45
     ...
     3.4938   57.8745    0.1888   56.9158;...  %iso 50
     2.8270   57.8910   -0.3643   58.0431;...  %iso 50
     3.1871   57.6251    0.1712   57.4410;...  %iso 50
     5.1890   58.0014   -1.1683   56.8878;...  %iso 50
     3.6037   57.4105    0.0698   58.5983;...  %iso 50
     ...
     4.2456   61.6943   -0.1278   60.1065;...  %iso 55
     4.1750   60.7749    0.0253   61.2117;...  %iso 55
     2.6828   60.7500   -1.0491   62.0551;...  %iso 55
     ...
     % 5.7915   67.3532   -0.4916   68.2567;...  %iso 65
     % 4.6038   68.8346    0.2200   69.6483;...  %iso 65
     % 4.9698   68.1026   -0.0030   68.4435;...  %iso 65
     ...
     4.9146   68.0136    0.0193   67.5259;...  %iso 65 (new data)
     4.1185   69.0463    3.3844   67.8286;...  %iso 65
     5.0331   67.7451    1.4734   68.2338;...  %iso 65
     5.2708   68.9428   -0.7598   68.3880;...  %iso 65
     4.4129   68.9013   -1.9455   69.3028;...  %iso 65 (new data)
     ...
     3.9979   72.6560    0.0242   73.8703;...  %iso 75
     5.3560   75.4868   -1.3309   74.0779;...  %iso 75
     3.5319   74.8857    3.4400   74.9433];    %iso 75
  
%% color pallete

ADilColor = [0.8500 0.3250 0.0980];
CDilColor = [0.9290 0.6940 0.1250];
ACColor = [0.4940 0.1840 0.5560];

SimColor = [0 0.4470 0.7410];

SimColorBase = [70 180 255]/255;

SimMap = [linspace(0,SimColorBase(1),18)', linspace(0,SimColorBase(2),18)', linspace(0,SimColorBase(3),18)'];
% figure; colormap(SimMap), colorbar;

SimColor75 = SimMap(end,:);
SimColor65 = SimMap(end-2,:);
SimColor55 = SimMap(end-4,:);
SimColor50 = SimMap(end-5,:);
SimColor45 = SimMap(end-6,:);
SimColor35 = SimMap(end-8,:);
SimColor0 = SimMap(end-15,:);

perpTextColor= [147 10 28]/255; parTextColor =[10 83 84]/255;

%% scatter the data plots
R_circ = 1000;
XC_circ = 0;
YC_circ = -R_circ;
x_circ_max = 20;
x_circ = linspace(-x_circ_max,x_circ_max,200);
y_circ = sqrt(R_circ^2-x_circ.^2)+YC_circ;

markerSize = 50; markeeralpha = 0.9;

Fig = figure('Position',[400, 200, 1150, 375]);
subplot(121); hold on, axis image;
fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');

s1 = scatter(data_pA(:,3), data_pA(:,1),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,3), data_pAC_dilute(:,1),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,3), data_pAC_dense(:,1),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');
% s4 = scatter(simDataResults(2:end,3), simDataResults(2:end,1),'filled','MarkerFaceColor',SimColor,'DisplayName','simulate');
% s4 = scatter(simDataGroupsResults(9:end,3), simDataGroupsResults(9:end,1),'filled','MarkerFaceColor',SimColor,'MarkerFaceAlpha',0.8,'DisplayName','simulate');

s50 = scatter(simDataGroupsResults(simDataGroups(:,1)==50,3), simDataGroupsResults(simDataGroups(:,1)==50,1),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','sim, \sigma=50');
s65 = scatter(simDataGroupsResults(simDataGroups(:,1)==65,3), simDataGroupsResults(simDataGroups(:,1)==65,1),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','sim, \sigma=65');


ax = gca;
ax.YAxisLocation= 'origin';
ax.XAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';


% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");

xlabel({"parallel", "mean"},'Position',[19,5,-1],'Color', parTextColor,'HorizontalAlignment','right');
ylabel("perpendicular mean",'Rotation',90,'Position',[1,-12.5,-1],'Color',perpTextColor);
xlim([-x_circ_max,x_circ_max]); ylim([-14,8]);
legend([s1,s2,s3,s50, s65],'Location',"southeast",'FontWeight','normal','FontSize',12,'LineWidth',0.5);

text(-19,-1.2,'dense','FontSize',12,'HorizontalAlignment','left');
text(-19,1.2,'dilute','FontSize',12,'HorizontalAlignment','left');


x_lin_min = 48;
x_lin_max = 82;
x_lin = linspace(x_lin_min,x_lin_max,100);
y_lin = x_lin;

subplot(122); hold on, axis image;
plot(x_lin,y_lin,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2);
s1 = scatter(data_pA(:,4), data_pA(:,2),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,4), data_pAC_dilute(:,2),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,4), data_pAC_dense(:,2),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');
% s4 = scatter(simDataResults(2:end,4), simDataResults(2:end,2),'filled','MarkerFaceColor',SimColor,'DisplayName','simulate');
% s4 = scatter(simDataGroupsResults(9:end,4), simDataGroupsResults(9:end,2),'filled','MarkerFaceColor',SimColor,'MarkerFaceAlpha',0.8,'DisplayName','simulate');
s50 = scatter(simDataGroupsResults(simDataGroups(:,1)==50,4), simDataGroupsResults(simDataGroups(:,1)==50,2),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','sim, \sigma=50');
s65 = scatter(simDataGroupsResults(simDataGroups(:,1)==65,4), simDataGroupsResults(simDataGroups(:,1)==65,2),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','sim, \sigma=65');


xlim([x_lin_min,x_lin_max]); ylim([x_lin_min,x_lin_max]);
xlabel("parallel \sigma",'Position',[56.7,51.2,-1],'Color', parTextColor);
ylabel("perpendicular \sigma",'Position',[51.2,60.5,-1],'Color', perpTextColor);
legend([s1,s2,s3,s50, s65],'Location','northwest','FontWeight','normal','FontSize',12,'LineWidth',0.5);

ax = gca;
ax.YAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';

%%
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures";
imgName =  "Scatter of the estimates";

exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);
%% scatter the data plots - center data
R_circ = 1000;
XC_circ = 0;
YC_circ = -R_circ;
x_circ_max = 20;
x_circ = linspace(-x_circ_max,x_circ_max,200);
y_circ = sqrt(R_circ^2-x_circ.^2)+YC_circ;

markerSize = 50; markeeralpha = 0.9;

Fig = figure('Position',[400, 200, 1150, 375]);
subplot(121); hold on, axis image;
fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');

s1 = scatter(data_pA(:,3), data_pA(:,1),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,3), data_pAC_dilute(:,1),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,3), data_pAC_dense(:,1),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');

s5 = scatter(data_pA_center(:,3), data_pA_center(:,1),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-interior');
s6 = scatter(data_pAC_center(:,3), data_pAC_center(:,1),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-interior');


ax = gca;
ax.YAxisLocation= 'origin';
ax.XAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';


% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");

xlabel({"parallel", "mean"},'Position',[19,5,-1],'Color', parTextColor,'HorizontalAlignment','right');
ylabel("perpendicular mean",'Rotation',90,'Position',[1,-12.5,-1],'Color',perpTextColor);
xlim([-x_circ_max,x_circ_max]); ylim([-14,8]);
legend([s1,s2,s3,s5, s6],'Location',"southeast",'FontWeight','normal','FontSize',12,'LineWidth',0.5);

text(-19,-1.2,'dense','FontSize',12,'HorizontalAlignment','left');
text(-19,1.2,'dilute','FontSize',12,'HorizontalAlignment','left');


x_lin_min = 48;
x_lin_max = 82;
x_lin = linspace(x_lin_min,x_lin_max,100);
y_lin = x_lin;

subplot(122); hold on, axis image;
plot(x_lin,y_lin,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2);
s1 = scatter(data_pA(:,4), data_pA(:,2),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,4), data_pAC_dilute(:,2),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,4), data_pAC_dense(:,2),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');
s5 = scatter(data_pA_center(:,4), data_pA_center(:,2),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-interior');
s6 = scatter(data_pAC_center(:,4), data_pAC_center(:,2),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-interior');



xlim([x_lin_min,x_lin_max]); ylim([x_lin_min,x_lin_max]);
xlabel("parallel \sigma",'Position',[56.7,51.2,-1],'Color', parTextColor);
ylabel("perpendicular \sigma",'Position',[51.2,60.5,-1],'Color', perpTextColor);
legend([s1,s2,s3,s5, s6],'Location','northwest','FontWeight','normal','FontSize',12,'LineWidth',0.5);

ax = gca;
ax.YAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';

%%
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures";
imgName =  "Scatter of the estimates - interface_interior";

exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);

%%
Fig = figure('Position',[400, 200, 400, 435]);

x_lin_min = -5;
x_lin_max = 95;
x_lin = linspace(x_lin_min,x_lin_max,100);
y_lin = x_lin;

hold on, axis image;
plot(x_lin,x_lin,"Color",[0.75 0.15 0.15],'LineStyle','--','LineWidth',2);

markerSize = 50; markeeralpha = 0.8;
% scatter(simData(1:end-1),(simDataResults(1:end-1,4)+simDataResults(1:end-1,2))/2,40,'filled','MarkerFaceColor',SimColor,'MarkerFaceAlpha',0.8,'DisplayName','sim');
s0 = scatter(simDataGroups(simDataGroups(:,1)==0),...
    (simDataGroupsResults(simDataGroups(:,1)==0,4)+simDataGroupsResults(simDataGroups(:,1)==0,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor0,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=0');
s35 = scatter(simDataGroups(simDataGroups(:,1)==35),...
    (simDataGroupsResults(simDataGroups(:,1)==35,4)+simDataGroupsResults(simDataGroups(:,1)==35,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=35');
s45 = scatter(simDataGroups(simDataGroups(:,1)==45),...
    (simDataGroupsResults(simDataGroups(:,1)==45,4)+simDataGroupsResults(simDataGroups(:,1)==45,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor45,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=45');
s50 = scatter(simDataGroups(simDataGroups(:,1)==50),...
    (simDataGroupsResults(simDataGroups(:,1)==50,4)+simDataGroupsResults(simDataGroups(:,1)==50,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor50,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=50');
s55 = scatter(simDataGroups(simDataGroups(:,1)==55),...
    (simDataGroupsResults(simDataGroups(:,1)==55,4)+simDataGroupsResults(simDataGroups(:,1)==55,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor55,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=55');
s65 = scatter(simDataGroups(simDataGroups(:,1)==65),...
    (simDataGroupsResults(simDataGroups(:,1)==65,4)+simDataGroupsResults(simDataGroups(:,1)==65,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor65,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=65');
s75 = scatter(simDataGroups(simDataGroups(:,1)==75),...
    (simDataGroupsResults(simDataGroups(:,1)==75,4)+simDataGroupsResults(simDataGroups(:,1)==75,2))/2,markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=75');

% newSimMap = [linspace(SimColor0(1), SimColor75(1), 16)', linspace(SimColor0(2), SimColor75(2), 16)',linspace(SimColor0(3), SimColor75(3), 16)'];
% hcb=colorbar('east','Position',[0.73,0.24,0.06,0.6]); colormap(newSimMap); clim([-2.5,77.5]);
% hcb.Title.String = 'GT \sigma';
% set(hcb,'XTick',0:5:75, 'XTickLabel',{'0','','','','','','','35','','45','50','55','','65','','75'},'fontWeight',"bold");

legend([s75, s65, s55, s50, s45, s35, s0],'Position',[0.64,0.27,0.23,0.32],'FontWeight','normal','FontSize',12,'LineWidth',0.5);

xlim([min(x_lin),max(x_lin)]); ylim([min(y_lin),max(y_lin)]);
xlabel("GT",'Position',[93,1,-1],'FontWeight','bold');
ylabel("Estimated",'Rotation',90,'Position',[2,66,-1],'FontWeight','bold');
set(gca,"XAxisLocation",'origin',"YAxisLocation",'origin','FontSize',12)
set(gca,"XTick",[0,35,55,75]);

ax = gca;
ax.YAxisLocation= 'origin';
ax.FontSize = 12;
ax.LineWidth = 1.2;
% ax.FontWeight = 'bold';

% legend();

%%
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures";
imgName =  "Simulate results - set displm and est displm";

exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);


%%

markerSize = 50; markeeralpha = 0.9;

Fig = figure('Position',[400, 200, 1150, 375]);
subplot(121); hold on, axis image;
fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');

s0 = scatter(simDataGroupsResults(simDataGroups(:,1)==0,3), simDataGroupsResults(simDataGroups(:,1)==0,1),markerSize,...
    'filled','MarkerFaceColor',SimColor0,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=0');
s35 = scatter(simDataGroupsResults(simDataGroups(:,1)==35,3), simDataGroupsResults(simDataGroups(:,1)==35,1),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=35');
s45 = scatter(simDataGroupsResults(simDataGroups(:,1)==45,3), simDataGroupsResults(simDataGroups(:,1)==45,1),markerSize,...
    'filled','MarkerFaceColor',SimColor45,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=45');
s55 = scatter(simDataGroupsResults(simDataGroups(:,1)==55,3), simDataGroupsResults(simDataGroups(:,1)==55,1),markerSize,...
    'filled','MarkerFaceColor',SimColor55,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=55');
s75 = scatter(simDataGroupsResults(simDataGroups(:,1)==75,3), simDataGroupsResults(simDataGroups(:,1)==75,1),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=75');
s50 = scatter(simDataGroupsResults(simDataGroups(:,1)==50,3), simDataGroupsResults(simDataGroups(:,1)==50,1),markerSize,...
    'filled','MarkerFaceColor',SimColor50,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=50');
s65 = scatter(simDataGroupsResults(simDataGroups(:,1)==65,3), simDataGroupsResults(simDataGroups(:,1)==65,1),markerSize,...
    'filled','MarkerFaceColor',SimColor65,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=65');


ax = gca;
ax.YAxisLocation= 'origin';
ax.XAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';


% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");

xlabel({"parallel", "mean"},'Position',[14.5,0,-1],'Color', parTextColor,'HorizontalAlignment','right');
ylabel("perpendicular mean",'Rotation',90,'Position',[0.8,-3.1,-1],'Color',perpTextColor);
xlim([-15,15]); ylim([-4,10]); yticks(-2:2:8);
legend([s75, s65, s55, s50, s45, s35, s0],'Location',"southoutside",'FontWeight','normal','FontSize',12,'LineWidth',0.5,'Orientation',"horizontal",'NumColumns',4);

text(-14.5,-0.8,'dense','FontSize',12,'HorizontalAlignment','left');
text(-14.5,0.8,'dilute','FontSize',12,'HorizontalAlignment','left');

% 
% xlabel("parallel mean",'Position',[13.5,-11,-1],'Color', parTextColor);
% ylabel("perpendicular mean",'Rotation',90,'Position',[1,-12.5,-1],'Color',perpTextColor);
% xlim([-x_circ_max,x_circ_max]); ylim([-14,8]);
% legend([s1,s2,s3,s4],'Location',"southwest",'FontWeight','normal','FontSize',12,'LineWidth',0.5);
% 
% text(19,-1.2,'dense','FontSize',12,'HorizontalAlignment','right');
% text(19,1.2,'dilute','FontSize',12,'HorizontalAlignment','right');



x_lin_min = -5;
x_lin_max = 85;
x_lin = linspace(x_lin_min,x_lin_max,100);
y_lin = x_lin;

subplot(122); hold on, axis image;
plot(x_lin,y_lin,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2);

s0 = scatter(simDataGroupsResults(simDataGroups(:,1)==0,4), simDataGroupsResults(simDataGroups(:,1)==0,2),markerSize,...
    'filled','MarkerFaceColor',SimColor0,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=0');
s35 = scatter(simDataGroupsResults(simDataGroups(:,1)==35,4), simDataGroupsResults(simDataGroups(:,1)==35,2),markerSize,...
    'filled','MarkerFaceColor',SimColor35,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=35');
s45 = scatter(simDataGroupsResults(simDataGroups(:,1)==45,4), simDataGroupsResults(simDataGroups(:,1)==45,2),markerSize,...
    'filled','MarkerFaceColor',SimColor45,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=45');
s55 = scatter(simDataGroupsResults(simDataGroups(:,1)==55,4), simDataGroupsResults(simDataGroups(:,1)==55,2),markerSize,...
    'filled','MarkerFaceColor',SimColor55,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=55');
s75 = scatter(simDataGroupsResults(simDataGroups(:,1)==75,4), simDataGroupsResults(simDataGroups(:,1)==75,2),markerSize,...
    'filled','MarkerFaceColor',SimColor75,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=75');

s50 = scatter(simDataGroupsResults(simDataGroups(:,1)==50,4), simDataGroupsResults(simDataGroups(:,1)==50,2),markerSize,...
    'filled','MarkerFaceColor',SimColor50,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=50');
s65 = scatter(simDataGroupsResults(simDataGroups(:,1)==65,4), simDataGroupsResults(simDataGroups(:,1)==65,2),markerSize,...
    'filled','MarkerFaceColor',SimColor65,'MarkerFaceAlpha',markeeralpha,'DisplayName','\sigma=65');


xlim([x_lin_min,x_lin_max]); ylim([x_lin_min,x_lin_max]);
xlabel("parallel \sigma",'Position',[34,2,-1],'Color', parTextColor);
ylabel("perpendicular \sigma",'Rotation',90,'Position',[2,13,-1],'Color', perpTextColor);
legend([s75, s65, s55, s50, s45, s35, s0],'Position',[0.78,0.20,0.076,0.37],'FontWeight','normal','FontSize',12,'LineWidth',0.5);

ax = gca;
ax.YAxisLocation= 'origin';
ax.XAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';

%%
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures";
imgName =  "Scatter of the estimates - only simulations";

exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);














%
%%
%% scatter the data plots - with all sim data above 45

markerSize = 50; markeeralpha = 0.9;

Fig = figure('Position',[400, 200, 1150, 375]);
subplot(121); hold on, axis image;
fill([x_circ, fliplr(x_circ)], [y_circ, -400*ones(size(y_circ))], [0.3 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');

s1 = scatter(data_pA(:,3), data_pA(:,1),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,3), data_pAC_dilute(:,1),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,3), data_pAC_dense(:,1),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');
% s4 = scatter(simDataResults(2:end,3), simDataResults(2:end,1),'filled','MarkerFaceColor',SimColor,'DisplayName','simulate');
s4 = scatter(simDataGroupsResults(9:end,3), simDataGroupsResults(9:end,1),markerSize,'filled','MarkerFaceColor',SimColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','simulate');


ax = gca;
ax.YAxisLocation= 'origin';
ax.XAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';


% plot(x_circ,y_circ,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2,'DisplayName',"Interface curvature");

xlabel({"parallel", "mean"},'Position',[19,5,-1],'Color', parTextColor,'HorizontalAlignment','right');
ylabel("perpendicular mean",'Rotation',90,'Position',[1,-12.5,-1],'Color',perpTextColor);
xlim([-x_circ_max,x_circ_max]); ylim([-14,8]);
legend([s1,s2,s3,s4],'Location',"southeast",'FontWeight','normal','FontSize',12,'LineWidth',0.5);

text(-19,-1.2,'dense','FontSize',12,'HorizontalAlignment','left');
text(-19,1.2,'dilute','FontSize',12,'HorizontalAlignment','left');

x_lin_min = 48;
x_lin_max = 82;
x_lin = linspace(x_lin_min,x_lin_max,100);
y_lin = x_lin;

subplot(122); hold on, axis image;
plot(x_lin,y_lin,'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth',2);
s1 = scatter(data_pA(:,4), data_pA(:,2),markerSize,...
    'filled','MarkerFaceColor',ADilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-dil');
s2 = scatter(data_pAC_dilute(:,4), data_pAC_dilute(:,2),markerSize,...
    'filled','MarkerFaceColor',CDilColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC-dil');
s3 = scatter(data_pAC_dense(:,4), data_pAC_dense(:,2),markerSize,...
    'filled','MarkerFaceColor',ACColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','pA-pC');
% s4 = scatter(simDataResults(2:end,4), simDataResults(2:end,2),'filled','MarkerFaceColor',SimColor,'DisplayName','simulate');
s4 = scatter(simDataGroupsResults(9:end,4), simDataGroupsResults(9:end,2),markerSize,'filled','MarkerFaceColor',SimColor,'MarkerFaceAlpha',markeeralpha,'DisplayName','simulate');


xlim([x_lin_min,x_lin_max]); ylim([x_lin_min,x_lin_max]);
xlabel("parallel \sigma",'Position',[56.7,51.2,-1],'Color', parTextColor);
ylabel("perpendicular \sigma",'Position',[51.2,60.5,-1],'Color', perpTextColor);
legend([s1,s2,s3,s4],'Location','northwest','FontWeight','normal','FontSize',12,'LineWidth',0.5);

ax = gca;
ax.YAxisLocation= 'origin';
ax.YColor = perpTextColor;
ax.XColor = parTextColor;
ax.FontSize = 12;
ax.LineWidth = 1.3;
ax.FontWeight = 'bold';

%%
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures";
imgName =  "Scatter of the estimates - with all sim data above 45";

exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);