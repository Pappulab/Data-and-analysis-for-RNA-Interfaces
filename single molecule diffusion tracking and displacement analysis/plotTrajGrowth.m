function plotTrajGrowth(loc_list,mol_select,mol_ind,dataNum,threshold,check_frames,Savepath,colors,preNaming,axisRang)
colorSz = size(colors,1);
Fig1 = figure('Position',[475,114,740,600]); hold on; axis image
for growthFactor = 1:15
    for i = 1:length(mol_select)
        mol_selectNum = mol_select(i);
        loc_select = loc_list(loc_list(:,mol_ind)==mol_selectNum,:);
        x_loc = loc_select(:,2);
        y_loc = loc_select(:,3);
        color_ind = mod(i,colorSz);
        color_ind(color_ind==0) = colorSz;
        colorSelect =  colors(color_ind,:);
        if length(x_loc)>growthFactor
            plot(x_loc(1:growthFactor+1), y_loc(1:growthFactor+1),'LineWidth',1.2,'Color',colorSelect); hold on;
            if growthFactor == 1
            scatter(x_loc(1,:), y_loc(1,:),15,colorSelect(1:3),'filled');
            end
        end
    end
    
    if axisRang~=[]
        xlim([axisRang(1), axisRang(2)]);
        ylim([axisRang(3), axisRang(4)]);
    else
        xlim([-1600, 1600]);
        ylim([-1600, 1600]);
    end
    saveresultDir = strcat(Savepath,"\traj growth\");
    if ~exist(saveresultDir, 'dir')
       mkdir(saveresultDir)
    end
    title(strcat("Trajectory - "," Growth = ", num2str(growthFactor)," ms"));
    % subtitle(strcat('Molecules [', num2str(mol_select),']'));
    xlabel('X position nm'); ylabel('Y position nm');
    set(gca,'FontSize',14)
    exportgraphics(Fig1,strcat(saveresultDir,"\threshold ", num2str(threshold)," checkFrames ", ...
                num2str(check_frames)," molecules ",preNaming," trajectory locs ", num2str(dataNum)," Growth = ", num2str(growthFactor)," ms", '.jpg'),'Resolution',600);

end
end