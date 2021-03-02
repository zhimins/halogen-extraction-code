clc
clear
warning off
foldername=uigetdir;
cd(foldername);
addpath(genpath(pwd));
list=dir([foldername,'\data\dataUMF\','*.csv']);%'*.csv'%% folder name for the UMF data stored
list_AMF=dir([foldername,'\data\dataAMF\','*.CSV']);%'*.csv'%% folder name for the AMF data stored

%% input the isotop want to check Cl or Br
element=input('input Cl or Br: ','s')
if element=='Cl'
    q=100/31.6  %% the itensity ratio for Cl35 and Cl37
    deltam=1.997 %% mass difference between Cl35 and Cl37
else
    q=0.9786  %% the itensity ratio for Br79 and Br81
    deltam=1.998  %% mass difference between Br79 and Br81
end   

% q=input('input the specific value of intensity that M/M+2, supposing that the intensity of M+2 is 1:')%%for cl :100/31.96; for Br :0.9786
% p=1;
% deltam=input('The difference in the mass of the isotope:') %% Cl is 1.997(Cl37-Cl35); Br is 1.998
%% input the ppm error allowed, 20 is used for song
percentage_error=input('the allowable error threshold for the specific value of intensity, the format is a decimal:') 

%%
[Cl1_0,Cl1_1,Cl2_0,ppm_allow1,ppm_allow2] = ppm_range(list,list_AMF,q,percentage_error,deltam,foldername);
[list_MS]=MS_A(foldername,Cl1_0,Cl1_1,Cl2_0,ppm_allow1,ppm_allow2,percentage_error,deltam);
[list_UMF]=UMF_unsureCl(list_MS,foldername);

