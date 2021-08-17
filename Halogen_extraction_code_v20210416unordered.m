clc
clear
warning off
foldername=uigetdir;
cd(foldername);
addpath(genpath(pwd));
list=dir([foldername,'\data\dataUMF\','*_UMF.CSV']);%'*_UMF.csv stored in '\foldername\data\dataUMF\'; dir means reading the files
list_AMF=dir([foldername,'\data\dataAMF\','*_AMF.CSV']);%'*_UMF.csv stored in '\foldername\data\dataAMF\'

%% input the mass difference and abundance ratio of isotop to check Cl or Br
Cl37_35=31.96/100;  %% the natural abundance ratio for Cl35 and Cl37
delta_m_Cl=1.997; %% mass difference between Cl35 and Cl37
Br81_79=0.9728;  %%  the natural abundance ratio for Br79 and Br81, source is MS textbook
delta_m_Br=1.998;  %% mass difference between Br79 and Br81
delta_m=[delta_m_Cl,delta_m_Br]; 

%% input the abundance error allowed, 0.2 is used for song,0.3 is used for liang
percentage_error=input('the allowable error, usually 0.3: '); 
recal_ppm=input('MFAssignR recalibration mass ppm, usually 10 or 13:');

%%
X0_1=[Cl37_35,Br81_79];%X0_1 is the abundance ratio for X(Cl37/Cl35 or Br81/Br79) in the sample
X1_1=[2*Cl37_35,2*Br81_79];% X1_1 is the abundance ratio for X(Cl35Cl37/Cl35Cl35 or Br79Br81/Br79Br79) in the sample
X2_1=[3*Cl37_35,3*Br81_79];% X2_1 is the abundance ratio for X(Cl35Cl35Cl37/Cl35Cl35Cl35 or Br79Br79Br81/Br79Br79Br79) in the sample
ppm_allow1=[20,20];% ppm_allow1 means X(Cl or Br)mass difference ppm error 
list_MS=dir([foldername,'\data\dataMS\align\blanksub\','*_MS.csv']);%'*_MS.csv stored in '\foldername\data\dataMS\'
h=waitbar(0,'MS_A is calculating');
for i=1:size(list_MS,1)%i=47,2021-04-16,9:34
    [list_MS]=MS_A_both(list_MS,i,foldername,X0_1(1),X1_1(1),X2_1(1),ppm_allow1(1),X0_1(2),...
        X1_1(2),X2_1(2),ppm_allow1(2),percentage_error,delta_m(1),delta_m(2),recal_ppm);%
    waitbar(i/size(list_MS,1))
end
close(h)
[list_UMF]=UMF_unsureX(list_MS,foldername);


