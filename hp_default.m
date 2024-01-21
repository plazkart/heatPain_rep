%This script loads mainStruct and sets main directory

%set path to folder in which all data is located
mainFolder = 'G:\_other\fMRS-heatPain';
load([mainFolder '\_meta\structMeta.mat']);

mainStruct.meta.folder = mainFolder;
spm_default_path = 'C:\Users\Science\Documents\MATLAB\spm12';
mainStruct.meta.SPMfolder = spm_default_path;