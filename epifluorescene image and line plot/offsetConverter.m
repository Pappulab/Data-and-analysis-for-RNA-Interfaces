tifName = 'C:\Users_NotBackedUp\Yuanxin\YQ - epifluoresence plots\data\offset_1\offset_cropped.tif';
tiff_info = imfinfo(tifName); % return tiff structure, one element per image
tiff_stack = double(imread(tifName, 1)) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = double(imread(tifName, ii));
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
offset = mean(tiff_stack,3);
%%
save('C:\Users_NotBackedUp\Yuanxin\YQ - epifluoresence plots\offSet.mat',"offset");