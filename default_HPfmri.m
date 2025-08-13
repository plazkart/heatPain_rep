%This script loads mainStruct and sets main directory

%set path to folder in which all data is located
mainFolder = 'E:\Alex\fmri_thermal';
load([mainFolder '\_meta\structMeta.mat']);

mainStruct.meta.folder = mainFolder;
spm_default_path = 'C:\Users\User\Documents\MATLAB';
mainStruct.meta.SPMfolder = spm_default_path;

% YD_default_path = 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results';
% mainStruct.meta.YDfolder = YD_default_path;