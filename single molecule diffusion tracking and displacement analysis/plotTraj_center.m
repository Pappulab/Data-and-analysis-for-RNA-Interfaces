function plotTraj_center(loc_list_sorted,mol_select,mol_ind,seq_jump_par_ind,seq_jump_perp_ind,dataNum,threshold,check_frames,path,colors,preNaming,axisRang)

colorSz = size(colors,1);
Fig1 = figure('Position',[475,114,740,600]); hold on; axis image
title(strcat("Trajectory centered - ",preNaming));
% subtitle(strcat('Molecules [', num2str(mol_select),']'));
xlabel('Parallel (nm)'); ylabel('Perpendicular (nm)');
% xlim([-(threshold*1.5+50),(threshold*1.5+50)]); ylim([-(threshold*1.5+50),(threshold*1.5+50)]);
set(gca,'FontSize',14)

for i = 1:length(mol_select)
    mol_selectNum = mol_select(i);
    loc_select = loc_list_sorted(loc_list_sorted(:,mol_ind)==mol_selectNum,:);
    par_jump = loc_select(:,seq_jump_par_ind);
    perp_jump = loc_select(:,seq_jump_perp_ind);
    par_loc = cumsum(par_jump);
    perp_loc = cumsum(perp_jump);
    color_ind = mod(i,colorSz);
    color_ind(color_ind==0) = colorSz;
    colorSelect =  colors(color_ind,:);
    plot(par_loc, perp_loc,'LineWidth',1.2,'Color',colorSelect); hold on;
end

if ~isempty(axisRang)
    xlim([axisRang(1), axisRang(2)]);
    ylim([axisRang(3), axisRang(4)]);
else
    xlim([-600, 600]);
    ylim([-600, 600]);
end
% exportgraphics(Fig1,strcat(path,"\threshold ", num2str(threshold)," checkFrames ", ...
%             num2str(check_frames)," molecules",preNaming," Centered trajectory locs ", num2str(dataNum),'.jpg'),'Resolution',600);
end