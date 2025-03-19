function [loc_list,mol_ind] = classifyDiffuse(loc_list,threshold,check_frames,frameN)
loc_list_dim = size(loc_list);
mol_ind = loc_list_dim(2)+1;
% multiMol_ind = loc_list_dim(2)+2;
loc_list = [loc_list, zeros(loc_list_dim(1),1)];%,zeros(loc_list_dim(1),1)];
molNum = 0;
for frameInd = 1:frameN-1
    frame1 = frameInd;
    loc_f1 = loc_list(loc_list(:,1)==frame1,2:3);
    % giving zero idx molecules index
    mol_f1 = loc_list(loc_list(:,1)==frame1,mol_ind);
%     multiMol_1 = loc_list(loc_list(:,1)==frame1,multiMol_ind);
    for idx = 1:length(mol_f1)
        if mol_f1(idx)==0
            mol_f1(idx)=molNum+1;
            molNum = molNum+1;
        else
        end
    end
    loc_list(loc_list(:,1)==frame1,mol_ind) = mol_f1;
%     loc_list(loc_list(:,1)==frame1,multiMol_ind) = multiMol_1;
     
    for frameCur = 1:min(check_frames-1,frameN-frameInd)   
        frame_cur = frameInd+frameCur;
        loc_cur = loc_list(loc_list(:,1)==frame_cur,2:3);
        mol_cur = loc_list(loc_list(:,1)==frame_cur,mol_ind);
%         multiMol_cur = loc_list(loc_list(:,1)==frame_cur,multiMol_ind);
        
        for iii = 1:size(loc_f1,1)
            if ismember(mol_f1(iii),mol_cur)
                continue
            else
            [CurInd, d] = crossDistance1D(loc_f1(iii,:),loc_cur,threshold);
            mol_cur(CurInd) = mol_f1(iii);
            if ~isempty(CurInd)
                loc_f1(iii,:) = (loc_f1(iii,:)*frameCur+loc_cur(CurInd,:))/(frameCur+1);
            end
            end
        end

        % [CurInd, f1Ind, d] = crossDistance(loc_f1,loc_cur,threshold);
        % mol_cur(CurInd) = mol_f1(f1Ind);
%         multiMol_1(f1Ind) = 1;
%         multiMol_cur(CurInd) = 1;
        loc_list(loc_list(:,1)==frame_cur,mol_ind) = mol_cur;
%         loc_list(loc_list(:,1)==frame1,multiMol_ind) = multiMol_1;
%         loc_list(loc_list(:,1)==frame_cur,multiMol_ind) = multiMol_cur;
    end
end

if loc_list(loc_list(:,1)==frameN,mol_ind)==0
    loc_list(loc_list(:,1)==frameN,mol_ind)=molNum+1;
end
end