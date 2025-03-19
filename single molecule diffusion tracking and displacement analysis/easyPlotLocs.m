close all;
 figure; scatter(SM_est_save_all(:,2),SM_est_save_all(:,3),3,'w','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); 
axis equal; hold on; axis off
molInd = 43;
scatter(loc_list_sorted(loc_list_sorted(:,5)==molInd,2),loc_list_sorted(loc_list_sorted(:,5)==molInd,3),20,'red','filled','LineWidth',15); 