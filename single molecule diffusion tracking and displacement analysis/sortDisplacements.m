function [loc_list_sorted, dist_ind, mol_angle_ind, seq_jump_dist_ind, ...
    seq_jump_angle_ind, seq_jump_par_ind,seq_jump_perp_ind, total_dist] = ...
    sortDisplacements(loc_list, mol_ind, XC, YC)
%%% Curent version 1: frame; 2-4: x, y, br; 5: mol ind; 6: traj length
%%% 7: mol interface angle
%%% 8: displacement; 
%%% 9: displacement azimuthal angle (average should be close to 0)
%%% 10-11: par/perp displacement
loc_list_sorted = sortrows(loc_list,mol_ind);
mol_index = unique(loc_list_sorted(:,mol_ind),'first');
mol_range = sort(mol_index,'ascend');
traj_length_frame = size(mol_range);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
dist_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
mol_angle_ind = size(loc_list_sorted,2);
loc_list_sorted(:,mol_angle_ind) = atan2d(loc_list_sorted(:,3)-YC,loc_list_sorted(:,2)-XC);
%(might work for evaluation but not used here or later yet)

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_dist_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_angle_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_par_ind = size(loc_list_sorted,2);

loc_list_sorted =  [loc_list_sorted,zeros(size(loc_list_sorted,1),1)];
seq_jump_perp_ind = size(loc_list_sorted,2);


total_dist=[];
for idx = 1:length(mol_range)
    molN = mol_range(idx);
    traj_length_frame(idx) = max(loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,1))-min(loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,1))+1;
    if traj_length_frame(idx)>1
        loc_temp =  sortrows(loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,:),1);
        %loc_temp =  loc_list_sorted(loc_list_sorted(:,mol_ind)==idx,:);
        % average jump (no angles)
        dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,dist_ind);
        
        for jj = 2:size(loc_temp,1)
            dist_temp = sqrt((loc_temp(jj-1,2)-loc_temp(jj,2))^2+(loc_temp(jj-1,3)-loc_temp(jj,3))^2);
            dist(jj) = dist(jj-1)+dist_temp;
           
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,dist_ind) = dist;
 
        total_dist = [total_dist;dist(end)];


        % sequential jumps (with angle)
        seq_jump_dist = loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_dist_ind);
        seq_jump_angle = loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_angle_ind);
        seq_jump_par = loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_par_ind);
        seq_jump_perp = loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_perp_ind);
        for jj = 2:size(loc_temp,1)
            vec_pos = [loc_temp(jj,2)-XC,loc_temp(jj,3)-YC]';
            vec_pos_norm = vec_pos/sqrt(sum(vec_pos.^2));
            vec_jump = [loc_temp(jj,2)-loc_temp(jj-1,2),loc_temp(jj,3)-loc_temp(jj-1,3)]';
            jump_temp = sqrt((loc_temp(jj-1,2)-loc_temp(jj,2))^2+(loc_temp(jj-1,3)-loc_temp(jj,3))^2);
            seq_jump_dist(jj) = jump_temp/(loc_temp(jj,1)-loc_temp(jj-1,1));
            seq_jump_angle(jj) = atan2d(vec_jump(2),vec_jump(1));
            seq_jump_perp(jj) = vec_pos_norm'*vec_jump;
            seq_jump_par(jj) = -[-vec_pos_norm(2),vec_pos_norm(1)]*vec_jump;
        end
        loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_dist_ind) = seq_jump_dist;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_angle_ind) = seq_jump_angle;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_par_ind) = seq_jump_par;
        loc_list_sorted(loc_list_sorted(:,mol_ind)==molN,seq_jump_perp_ind) = seq_jump_perp;
    end
end
end
