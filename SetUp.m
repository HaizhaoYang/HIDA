function SetUp()

file_path = mfilename('fullpath');
tmp = strfind(file_path,'SetUp');
file_path = file_path(1:(tmp(end)-1));

addpath(genpath([file_path 'Diffusion']));


end
