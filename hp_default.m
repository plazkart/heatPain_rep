%This script loads mainStruct and sets main directory

%set path to folder in which all data is located
mainFolder = '\\WIN-JC5203G34U7\fMRS-heatPain';
load([mainFolder '\_meta\structMeta.mat']);

mainStruct.meta.folder = mainFolder;
spm_default_path = 'C:\Users\Science\Documents\MATLAB\spm12';
mainStruct.meta.SPMfolder = spm_default_path;

YD_default_path = 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results';
mainStruct.meta.YDfolder = YD_default_path;