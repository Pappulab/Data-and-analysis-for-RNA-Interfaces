function Fig1 = plotTraj_angle(loc_list_sorted,mol_select,mol_ind,seq_jump_par_ind,seq_jump_perp_ind,perpColor,parColor,preNaming,axisRange)
cMap = [linspace(parColor(1),perpColor(1),255)',...
        linspace(parColor(2),perpColor(2),255)',...
        linspace(parColor(3),perpColor(3),255)'];
Fig1 = figure('Position',[207.8571428571428,83.28571428571428,1080,755]); hold on; axis image
title(strcat("Trajectory: ",preNaming));
xlabel('X position nm'); ylabel('Y position nm');
set(gca,'FontSize',14)
theta_dic = linspace(0, 90, size(cMap,1));

cMapBar = [cMap, 0.7*ones(size(cMap,1),1)];

for i = 1:length(mol_select)
    mol_selectNum = mol_select(i);
    loc_select = loc_list_sorted(loc_list_sorted(:,mol_ind)==mol_selectNum,:);
    x_loc = loc_select(:,2);
    y_loc = loc_select(:,3);
    parJump = loc_select(:,seq_jump_par_ind);
    perpJump = loc_select(:,seq_jump_perp_ind);
    theta_list = atan2d(perpJump,parJump);
    theta_list(theta_list<=0&theta_list>=-180)=-theta_list(theta_list<=0&theta_list>=-180);
    theta_list(theta_list<=180&theta_list>=90)=180-theta_list(theta_list<=180&theta_list>=90);
    
    for ind = 2: size(x_loc,1)
        [~, color_ind] = min(abs(theta_list(ind)-theta_dic));
        plot(x_loc(ind-1:ind), y_loc(ind-1:ind),'LineWidth',2,'Color',cMapBar(color_ind,:)); hold on;
        if ind == 2
            scatter(x_loc(1,:), y_loc(1,:),50,'w','filled','MarkerFaceAlpha',0.9);
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

plot([min(0.72*axisRange(2),axisRange(2)-600),min(0.72*axisRange(2),axisRange(2)-600)+500], [axisRange(3)+100, axisRange(3)+100],'w','LineWidth', 4);
text(min(0.72*axisRange(2),axisRange(2)-600)+250,axisRange(3)+160,'500 nm','Color','w','FontSize',16,'FontWeight','bold','HorizontalAlignment', 'center');
set(gca,'Color','k');
set(gca,'XTickLabel',[],'YTickLabel',[]);

colormap(cMap);
c = colorbar;
c.Ticks = [0 1];
c.TickLabels = ["parallel", "perpendicular"];
c.Ruler.TickLabelRotation=25;

end
