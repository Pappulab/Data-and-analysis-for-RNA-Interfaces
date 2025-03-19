function Fig1 = plotTraj(loc_list,mol_select,mol_ind,colors,preNaming,axisRange)
colorSz = size(colors,1);
Fig1 = figure('Position',[207.8571428571428,83.28571428571428,1080,755]); hold on; axis image
title(strcat("Trajectory: ",preNaming));
xlabel('X position nm'); ylabel('Y position nm');
set(gca,'FontSize',14)

for i = 1:length(mol_select)
    mol_selectNum = mol_select(i);
    loc_select = loc_list(loc_list(:,mol_ind)==mol_selectNum,:);
    x_loc = loc_select(:,2);
    y_loc = loc_select(:,3);
    color_ind = mod(i,colorSz);
    color_ind(color_ind==0) = colorSz;
    colorSelect =  colors(color_ind,:);
    plot(x_loc, y_loc,'LineWidth',2,'Color',colorSelect); hold on;
    scatter(x_loc(1,:), y_loc(1,:),50,colorSelect(1:3),'filled','MarkerFaceAlpha',0.9);
end

if ~isempty(axisRange)
    axisRange = axisRange;
else
    axisRange = [-1600 1600 -1600 1600];
end
xlim([axisRange(1), axisRange(2)]);
ylim([axisRange(3), axisRange(4)]);

plot([min(0.72*axisRange(2),axisRange(2)-600),min(0.72*axisRange(2),axisRange(2)-600)+500], [axisRange(3)+100, axisRange(3)+100],'w','LineWidth', 4);
text(min(0.72*axisRange(2),axisRange(2)-600)+250,axisRange(3)+160,'500 nm','Color','w','FontSize',16,'FontWeight','bold','HorizontalAlignment', 'center');
set(gca,'Color','k');
set(gca,'XTickLabel',[],'YTickLabel',[]);

end
