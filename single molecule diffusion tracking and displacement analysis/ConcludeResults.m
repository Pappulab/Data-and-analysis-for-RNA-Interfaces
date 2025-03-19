clear;
close all;

MatResultsSavePath = "C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\mat_save\";
ImgResultsSavePath = "C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\fig_save\";

dataNum = 8;
FoVTotal = 6;

locNum = zeros(FoVTotal,1);
molNum = zeros(FoVTotal,1);
molNum_dil = zeros(FoVTotal,1);
molNum_den = zeros(FoVTotal,1);

jumpNum_dil = zeros(FoVTotal,1);
jumpNum_den = zeros(FoVTotal,1);

fitJump_dil = zeros(FoVTotal,4);
fitJump_den = zeros(FoVTotal,4);
br_all = [];
meanBr = zeros(FoVTotal,1);

% burstTime_dil = [];
% burstTime_den = []; % not yet implement

for i = 1:FoVTotal
    fovNum=i;

    MatFileDir = strcat(MatResultsSavePath,'data',num2str(dataNum),'_FoV',num2str(fovNum),'\');

    allMat = strcat(MatFileDir,'data',num2str(dataNum),'_FoV',num2str(fovNum),"_all.mat");
    DenMat = strcat(MatFileDir,'data',num2str(dataNum),'_FoV',num2str(fovNum),"_Dense_Dense_Interface_intermediate.mat");
    DilMat = strcat(MatFileDir,'data',num2str(dataNum),'_FoV',num2str(fovNum),"_Dil_Dense_Interface_intermediate.mat");
    
    load(allMat,...
        "loc_list",'molecule_numel','molecule_numel_all');


    locNum(i) = size(loc_list,1);
    molNum(i) = size(molecule_numel_all,1);

    molNum_den(i) = size(molecule_numel,1);

    molNum_dil(i) = molNum(i) - molNum_den(i);

    br_curr = molecule_numel_all(:,4);

    br_all = [br_all; br_curr];

    meanBr(i) = mean(br_curr);
    
    load(DilMat,'jump_list','pd_par','pd_perp');

    jumpNum_dil(i) = size(jump_list,1);
    fitJump_dil(i,:) = [pd_perp.mu, pd_perp.std, pd_par.mu, pd_par.std];

    load(DenMat,'jump_list','pd_par','pd_perp');

    jumpNum_den(i) = size(jump_list,1);
    fitJump_den(i,:) = [pd_perp.mu, pd_perp.std, pd_par.mu, pd_par.std]; 

    

end
meanBr_all = mean(br_all);

MatFileName = strcat('data',num2str(dataNum),'_allCombined');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));

%%
clear;
close all;

MatResultsSavePath = "C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\mat_save\";
ImgResultsSavePath = "C:\Users_NotBackedUp\Yuanxin\YQ - Localization data diffusion analysis\fig_save\";

dataNum = 2;
FoVTotal = 6;

locNum = zeros(FoVTotal,1);
molNum = zeros(FoVTotal,1);
jumpNum = zeros(FoVTotal,1);

fitJump = zeros(FoVTotal,4);
br_all = [];
meanBr = zeros(FoVTotal,1);

% burstTime_dil = [];
% burstTime_den = []; % not yet implement

for i = 1:FoVTotal
    fovNum=i;

    MatFileDir = strcat(MatResultsSavePath,'data',num2str(dataNum),'_FoV',num2str(fovNum),'\');

    allMat = strcat(MatFileDir,'data',num2str(dataNum),'_FoV',num2str(fovNum),"_all.mat");
    
    load(allMat,...
        "loc_list",'molecule_numel_all',"jump_list",'pd_par','pd_perp');


    locNum(i) = size(loc_list,1);
    molNum(i) = size(molecule_numel_all,1);

    br_curr = molecule_numel_all(:,4);

    br_all = [br_all; br_curr];

    meanBr(i) = mean(br_curr);

    jumpNum(i) = size(jump_list,1);
    fitJump(i,:) = [pd_perp.mu, pd_perp.std, pd_par.mu, pd_par.std];



end
meanBr_all = mean(br_all);

MatFileName = strcat('data',num2str(dataNum),'_allCombined');
save(strcat(MatResultsSavePath,MatFileName,'.mat'));