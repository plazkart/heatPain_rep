%This script loads mainStruct and sets main directory

%set path to folder in which all data is located
mainFolder = 'E:\_other\fMRS-heatPain';
load([mainFolder '\_meta\structMeta.mat']);

mainStruct.meta.folder = mainFolder;