function plotTrajGrowth_colorCoded(loc_list,mol_select,mol_ind,dataNum,threshold,check_frames,Savepath,GrowFrames,preNaming,axisRange)
Fig1 = figure('Position',[475,114,740,600]); hold on; axis image
init_color = [0, 0, 0];
last_color = [1, 0, 0];
colors = [linspace(init_color(1),last_color(1),GrowFrames)', linspace(init_color(2),last_color(2),GrowFrames)', linspace(init_color(3),last_color(3),GrowFrames)', 0.6*ones(GrowFrames,1)];
for growthFactor = 1:GrowFrames
    for i = 1:length(mol_select)
        mol_selectNum = mol_select(i);
        loc_select = loc_list(loc_list(:,mol_ind)==mol_selectNum,:);
        x_loc = loc_select(:,2);
        y_loc = loc_select(:,3);
        color_ind = growthFactor;
        colorSelect = colors(color_ind,:);
        if length(x_loc)>growthFactor
            plot(x_loc(growthFactor:growthFactor+1), y_loc(growthFactor:growthFactor+1),'LineWidth',1.2,'Color',colorSelect); hold on;
            if growthFactor == 1
            scatter(x_loc(1,:), y_loc(1,:),15,colorSelect(1:3),'filled');
            end
        end
    end
    
    if ~isempty(axisRange)
    axisRange = axisRange;
    else
        axisRange = [-1600 1600 -1600 1600];
    end
    xlim([axisRange(1), axisRange(2)]);
    ylim([axisRange(3), axisRange(4)]);

    saveresultDir = strcat(Savepath,"traj growth");
    if ~exist(saveresultDir, 'dir')
       mkdir(saveresultDir)
    end
    title(strcat("T = ", num2str(growthFactor*10)," ms"));
    % subtitle(strcat('Molecules [', num2str(mol_select),']'));
    % xlabel('X position nm'); ylabel('Y position nm');
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(gca, 'TickLength',[0 0])
    set(gca,'FontSize',14)
    exportgraphics(Fig1,strcat(saveresultDir,"\threshold ", num2str(threshold)," checkFrames ", ...
                num2str(check_frames)," molecules ",preNaming," trajectory locs ", num2str(dataNum)," Growth = ", num2str(growthFactor)," ms", '.jpg'),'Resolution',600);

end