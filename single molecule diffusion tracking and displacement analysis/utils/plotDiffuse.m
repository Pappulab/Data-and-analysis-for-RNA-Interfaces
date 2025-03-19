function plotDiffuse(loc_list,mol_select,mol_ind,cropImg,dataNum,threshold,check_frames)
loc_select = loc_list(loc_list(:,mol_ind)==mol_select,:);
frame_first = min(loc_select(:,1));
frame_last = max(loc_select(:,1));
loc_select(:,2:3)=loc_select(:,2:3)/58.5+45;
figure('Position',[100,100,700,300]);
I_sz = size(cropImg,1);

for frameNum = frame_first:frame_last
    SM_cur_loc = [];
    imagesc(cropImg(:,:,frameNum)); axis image; axis off; colorbar; hold on; clim([95,145])
    idx = find(loc_select(:,1)==frameNum);
    if ~isempty(idx)
        SM_cur_loc = loc_select(idx,2:3);
        scatter(SM_cur_loc(1),SM_cur_loc(2),50,'rx','LineWidth',2);
        scatter(SM_cur_loc(1)+I_sz,SM_cur_loc(2),50,'rx','LineWidth',2);
        line([I_sz+0.5,I_sz+0.5],[0,I_sz],'LineWidth',2,'color','w');
        if frameNum==frame_first
            loc_temp = SM_cur_loc;
        elseif frameNum>=frame_first+1
            loc_temp = [loc_temp;SM_cur_loc]; 
        end

%scatter(0+45,0+45,50,'rx','LineWidth',2);
    else
    end
    if frameNum>=frame_first+1
     line(loc_temp(:,1) ,loc_temp(:,2),'Color','k','LineWidth',2);  
    end
    frame1 = getframe(gcf);
    im1 = frame2im(frame1);
    [A1,map1] = rgb2ind(im1,256);
    if frameNum == frame_first
        imwrite(im1,['threshold ', num2str(threshold),' checkFrames ', num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff')
    else 
        imwrite(im1,['threshold ', num2str(threshold),' checkFrames ', num2str(check_frames),' molecule ',num2str(mol_select),' trajactory raw and locs ', num2str(dataNum),'.tif'],'tiff','WriteMode','append');
    end
end
end