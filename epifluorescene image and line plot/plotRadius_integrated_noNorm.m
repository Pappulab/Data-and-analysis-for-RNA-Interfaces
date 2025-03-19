clc;
clear;

load('offSet.mat');

Dig2Ph=0.29;

ImageDataPath = "C:\Users_NotBackedUp\Yuanxin\YQ - epifluoresence plots\data";
figureSaveDirectory = "C:\Users\q.yuanxin\Box\LewLab_Yuanxin\Projects\RNA condensates - poly-rA interface diffusion\figures\panelB\updated02242025";

% dataNum = 2:12;
% laserPower =  [0.1, 1, 10, 25, 50, 80, 50, 25, 10, 1, 0.1];
% exposureTime = 20; 
% Radius_set = 49;

% dataNum = 17:31;
% laserPower = [0.1, 1, 2, 5, 10, 25, 50, 80, 50, 25, 10, 5, 2, 1, 0.1];

% dataNum = [18,21];
% laserPower = [1,10];
% exposureTime = 10; 
% Radius_set= 30;

dataNum = 23;
laserPower = 50;
exposureTime = 10; 
Radius_set= 30;
% 
% dataNum = 18;
% laserPower = 1;
% exposureTime = 10; 
% Radius_set= 30;

% dataNum = 21;
% laserPower = 10;
% exposureTime = 10; 
% Radius_set= 30;

% XcenterLoc = [1576,500]; % 1-14
% % YcenterLoc = [471,502];
% fovSize = 181;

XcenterLoc = [1494,419]; % 16-33
fovSize = 101;

resizefactor = 1;

%%
for i=1:length(dataNum)
    if laserPower(i)<45 %%
        %%
        tifName = strcat(ImageDataPath,"\_",num2str(dataNum(i)),'\_',num2str(dataNum(i)),'_MMStack_Default.ome.tif');
        tiff_info = imfinfo(tifName); % return tiff structure, one element per image
        tiff_stack = double(imread(tifName, 1)) ; % read in first image
        %concatenate each successive tiff to tiff_stack
        for ii = 2 : size(tiff_info, 1)
            temp_tiff = double(imread(tifName, ii));
            tiff_stack = cat(3 , tiff_stack, temp_tiff);
        end
        dataTif = double(tiff_stack);
        avgDataTif = mean(dataTif,3)-offset;
        avgDataTif = avgDataTif*Dig2Ph;
        %% Cropping and computing
        cropData = avgDataTif(XcenterLoc(2)+1-(fovSize-1)/2: XcenterLoc(2)+1+(fovSize-1)/2,...
            XcenterLoc(1)+1-(fovSize-1)/2:XcenterLoc(1)+1+(fovSize-1)/2);
        centerLoc_new = (fovSize+1)/2;
        [x,y]=meshgrid(-(fovSize-1)/2:(fovSize-1)/2);
        radiusData = sqrt(x.^2+y.^2);
        radiusNorm = radiusData.*58.5;
        Radius_sample = (0:0.05:2700/58.5/30)*Radius_set*58.5;
        radius_ind =  (Radius_sample(2:end)+Radius_sample(1:end-1))/2;
        I_mean = zeros(size(radius_ind));
        I_std =  zeros(size(radius_ind));
        I_med =  zeros(size(radius_ind));
        I_90 = zeros(size(radius_ind));
        I_10 =  zeros(size(radius_ind));
        for ii = 1:length(Radius_sample)-1
            DataPc = cropData(radiusNorm>=Radius_sample(ii) & radiusNorm<Radius_sample(ii+1));
            I_mean(ii) = mean(DataPc,'all');
            I_std(ii) =  std(DataPc,0,'all');
            I_med(ii) = prctile(DataPc,50,'all');
            I_90(ii) = prctile(DataPc,90,'all');
            I_10(ii) =  prctile(DataPc,10,'all');
        end
        %%
        Fig = figure('Position',[600,300,323,424]); tiledlayout('flow'); hold on;
        subplot(211,'Position',[0.07,0.456,0.95,0.5]); 
        imagesc(imresize(cropData,resizefactor,'method','nearest'));  colormap gray; hcb = colorbar; axis image; box off; axis off; caxis([0,175]);
        line(resizefactor*[centerLoc_new centerLoc_new+1.5*Radius_set+5],resizefactor*[centerLoc_new centerLoc_new],'Color', 'w','LineStyle','--','LineWidth',1);
        if i == 1
        line(resizefactor*[fovSize-20 fovSize-20+17.0940],resizefactor*[fovSize-10 fovSize-10],'Color','w','LineWidth',2);
        text(resizefactor*(fovSize-25),resizefactor*(fovSize-20),'1\mum','Color','w','FontWeight','bold');
        end
        
        set(get(hcb,'Title'),'String','photons','fontsize',10, 'Color','k'); 
        
        hcb.TickLabelInterpreter='latex';
        hcb.FontSize = 11;
        hcb.Color = 'k';

        subplot(212,'Position',[0.18,0.11,0.68,0.32]); hold on; xlabel('Radius (in nm)','Interpreter','latex'); ylabel('photons');
        
        plot(radius_ind, I_mean,'DisplayName', 'Average','LineWidth',2,'Color',[0 0.4470 0.7410]); 
        if  max(I_mean+I_std)>220
            ylim([0,min(20*(ceil(max(I_mean+I_std)./20+1)))]);
        else 
            ylim([0,200]);
        end
        xlim([0,2700]); set(gca,'XTick',[0,500,1000,1500,2000,2500],'TickLength',[0.04 0.05]);
        fill([radius_ind, fliplr(radius_ind)], [I_mean-I_std, fliplr(I_mean+I_std)], ...
            [0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.3);
        ax = gca;
        set(ax,'FontSize',11,'Color','none','XColor','k', 'YColor','k');       
        ax.TickLabelInterpreter='latex';
        ax.FontSize = 11;
%%
        imgName = strcat("EpiFluoRadiusPlot_data",num2str(dataNum(i)),"_Excitation_",num2str(laserPower(i)),"_expose_",num2str(exposureTime),"ms");

        savefig(Fig,strcat(figureSaveDirectory,'\',imgName,'.fig'));
        exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);

    else
        %%
        table_name = strcat(ImageDataPath,"\data",num2str(dataNum(i)),'_xch.csv');
        T=readtable(table_name);
        
        loc_list = T{:,1:5};
        img_sz = fovSize; pix_sz = 1; bin_sz = 1;
        grid = 0 : bin_sz/pix_sz : (img_sz) *pix_sz;
        [X,Y] = meshgrid(grid,grid);
        I_loc = zeros(size(X));
        for ii = 1:length(grid)
            for j = 1:length(grid)
                locs_temp = [];
                locs_temp = find( loc_list(:,2)>=X(ii,j)-bin_sz/2 & ...
                                  loc_list(:,2)<X(ii,j)+bin_sz/2 &...
                                  loc_list(:,3)>=Y(ii,j)-bin_sz/2 &...
                                  loc_list(:,3)<Y(ii,j)+bin_sz/2);
                I_loc(ii,j) = numel(locs_temp);
            end
        end
        x = loc_list(:,2);
        y = loc_list(:,3);
        xCenter = (img_sz+1)/2;
        yCenter = (img_sz+1)/2;
        radiusData = sqrt((x-xCenter).^2+(y-yCenter).^2);
        radiusNorm = radiusData.*58.5;    
        Radius_sample = (0:0.008:2700/58.5/30)*Radius_set*58.5;
        radius_ind =  (Radius_sample(2:end)+Radius_sample(1:end-1))/2;
        LocNum = zeros(size(radius_ind));
        for ii = 1:length(Radius_sample)-1
            DataPc = find(radiusNorm>=Radius_sample(ii) & radiusNorm<Radius_sample(ii+1));
            LocNum(ii) = numel(DataPc);
            LocNum(ii) = numel(DataPc);%/(2*2*pi*radius_ind(ii))
        end

        % FWHM
        HmaxNum = max(LocNum)/2;

        leftIndexR = find(LocNum >= HmaxNum, 1, 'first');
        leftRad = radius_ind(leftIndexR-1) +  (HmaxNum-LocNum(leftIndexR-1)) * (radius_ind(leftIndexR)-radius_ind(leftIndexR-1))/(LocNum(leftIndexR)-LocNum(leftIndexR-1));


        rightIndexL = find(LocNum >= HmaxNum, 1, 'last');
        rightRad = radius_ind(rightIndexL) +  (HmaxNum-LocNum(rightIndexL)) * (radius_ind(rightIndexL+1)-radius_ind(rightIndexL))/(LocNum(rightIndexL+1)-LocNum(rightIndexL));
        
        FWHM = rightRad-leftRad;


        Fig = figure('Position',[600,300,323,424]); tiledlayout('flow'); hold on; 
        subplot(211,'Position',[0.07,0.456,0.95,0.5]); imagesc(imresize(I_loc,resizefactor,'method','nearest')); 
        colormap gray; hcb = colorbar; axis image; box off; axis off;  caxis([0,5]);
        line(pix_sz/bin_sz*resizefactor*[xCenter xCenter+1.5*Radius_set+5],pix_sz/bin_sz*resizefactor*[yCenter yCenter],'Color', 'w','LineStyle','--','LineWidth',1);

        % line([60 60+17.0940],[75 75],'Color','w','LineWidth',4);
        % text(61,79,'1um','Color','w','FontWeight','bold');
    
        set(get(hcb,'Title'),'String','locs','fontsize',10, 'Color','k'); 
        hcb.TickLabelInterpreter='latex';
        hcb.FontSize = 11;
        hcb.Color = 'k';
        subplot(212,'Position',[0.18,0.11,0.68,0.32]); hold on; xlabel('Radius (in nm)','Interpreter','latex'), ylabel('Summed Locs');
        plot(radius_ind, LocNum,'DisplayName', 'Average','LineWidth',2,'Color',[0 0.4470 0.7410]); ylim([0,160]); 
        xlim([0,2700]); set(gca,'XTick',[0,500,1000,1500,2000,2500],'TickLength',[0.04 0.05]);
        plot([leftRad, leftRad],[0, HmaxNum*2],'--','LineWidth',1.5,'Color',[0.8 0.2 0.2]);
        plot([rightRad, rightRad],[0, HmaxNum*2],'--','LineWidth',1.5,'Color',[0.8 0.2 0.2]);
        plot([leftRad, rightRad],[HmaxNum, HmaxNum],'-','LineWidth',1.5,'Color',[0.8 0.2 0.2]);
        text(leftRad-650, HmaxNum+10,{"FWHM",num2str(FWHM)+" nm"},'Color',[0.8 0.2 0.2],'HorizontalAlignment','center');
        set(gca,'FontSize',11,'Color','none','XColor','k', 'YColor','k');
        

        imgName = strcat("EpiFluoRadiusLocsPlot_data",num2str(dataNum(i)),"_Excitation_",num2str(laserPower(i)),"_expose_",num2str(exposureTime),"ms");

        savefig(Fig,strcat(figureSaveDirectory,'\',imgName,'.fig'));
        exportgraphics(Fig,strcat(figureSaveDirectory,'\',imgName,'.pdf'),'ContentType','vector','Resolution',600);

    end
end
