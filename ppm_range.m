function [Cl1_0,Cl1_1,Cl2_0,ppm_allow1,ppm_allow2] = ppm_range(list,list_AMF,q,percentage_error,deltam,foldername)
%C11_0 is the result of intensity Cl37/Cl35 with 1 Cl;  C11_1 is the result
%of intensity Cl37/Cl35 with 2 Cl;C12_0 is the result of intensity
%Cl37/Cl35 with 3 Cl;ppm_allow1 is the lower bonds of m/z;ppm_allow2 is the upper bonds of m/z

cal_numberj1=0; 
cal_numberj2=0; 
cal_numberj3=0; 
UMF_cal1=[];
UMF_cal2=[];
UMF_cal3=[];
cal_number=0;
for i=1:length(list)
    [~,~,raw]=xlsread([foldername,'\data\dataUMF\',list(i,1).name]);
    UMF_select=[];
    UMF_select1=[];
    if(size(list_AMF,1)>1)
        splitname_AMF = strtok(list(i,1).name,'.');
        splitname_AMF_=[splitname_AMF(1:size(splitname_AMF,2)-3),'AMF.csv'];
        for j=1:length(list_AMF)
            if(list_AMF(j,1).name(1:min([size(list_AMF(j,1).name,2),size(splitname_AMF_,2)]))==splitname_AMF_(1:min([size(list_AMF(j,1).name,2),size(splitname_AMF_,2)])))
                [~,~,raw_AMF]=xlsread([splitname_AMF(1:size(splitname_AMF,2)-3),'AMF.csv']);
                AMF_select=[];
                for k=2:size(raw_AMF,1)
                    if(0<cell2mat(raw_AMF(k,16))&&cell2mat(raw_AMF(k,18))==0)
                        AMF_select=[AMF_select;raw_AMF(k,:)];
                    end
                end
                UMF_select=[UMF_select;AMF_select];
                break
            end
        end
    end
    for j=2:size(raw,1)
        if(0<cell2mat(raw(j,16))&&cell2mat(raw(j,18))==0)
            UMF_select=[UMF_select;raw(j,:)];
        elseif(0<cell2mat(raw(j,16))+cell2mat(raw(j,18)))
            UMF_select1=[UMF_select1;raw(j,:)];
        end
    end
    A=ones(size(UMF_select,1),size(UMF_select1,1));
    for j=1:size(UMF_select,1)
        cal_numberk=0;
        for k=1:size(UMF_select1,1)
            A(j,k)=sum(abs(cell2mat(UMF_select(j,[6:11,20]))-cell2mat(UMF_select1(k,[6:11,20]))));
            if(A(j,k)<0.000001&&cell2mat(UMF_select(j,16))+cell2mat(UMF_select(j,18))==cell2mat(UMF_select1(k,16))+cell2mat(UMF_select1(k,18)))
                cal_numberk=cal_numberk+1;
                if(UMF_select1{k,18}==1&&UMF_select1{k,16}==0)
                    cal_numberj1=cal_numberj1+1;
                    UMF_cal1(1,cal_numberj1)=UMF_select1{k,1}/UMF_select{j,1};
                elseif(UMF_select1{k,18}==1&&UMF_select1{k,16}==1)
                    cal_numberj2=cal_numberj2+1;
                    UMF_cal2(1,cal_numberj2)=UMF_select1{k,1}/UMF_select{j,1};
                elseif(UMF_select1{k,18}==1&&UMF_select1{k,16}==2)
                    cal_numberj3=cal_numberj3+1;
                    UMF_cal3(1,cal_numberj3)=UMF_select1{k,1}/UMF_select{j,1};
                end
            end
        end
    end
end
UMF_cal1=[UMF_cal1,1/q];
UMF_cal2=[UMF_cal2,2/q];
UMF_cal3=[UMF_cal3,3/q];
UMF_cal1=sort(UMF_cal1,2);
UMF_cal2=sort(UMF_cal2,2);
UMF_cal3=sort(UMF_cal3,2);
Cl1_0=mean(UMF_cal1(max(round(length(UMF_cal1)*0.05),1):round(length(UMF_cal1)*0.95)));
Cl1_1=mean(UMF_cal2(max(round(length(UMF_cal2)*0.05),1):round(length(UMF_cal2)*0.95)));
Cl2_0=mean(UMF_cal3(max(round(length(UMF_cal3)*0.05),1):round(length(UMF_cal3)*0.95)));

for i=1:length(list)
    [~,~,raw]=xlsread([foldername,'\data\dataUMF\',list(i,1).name]);    
    UMF_select=[];
    UMF_select1=[];
    if(size(list_AMF,1)>1)
        splitname_AMF = strtok(list(i,1).name,'.');
        splitname_AMF_=[splitname_AMF(1:size(splitname_AMF,2)-3),'AMF.csv'];
        for j=1:length(list_AMF)
            if(list_AMF(j,1).name(1:min([size(list_AMF(j,1).name,2),size(splitname_AMF_,2)]))==splitname_AMF_(1:min([size(list_AMF(j,1).name,2),size(splitname_AMF_,2)])))
                [~,~,raw_AMF]=xlsread([splitname_AMF(1:size(splitname_AMF,2)-3),'AMF.csv']);
                AMF_select=[];
                for k=2:size(raw_AMF,1)
                    if(0<cell2mat(raw_AMF(k,16))&&cell2mat(raw_AMF(k,18))==0)
                        AMF_select=[AMF_select;raw_AMF(k,:)];
                    end
                end
                UMF_select=[UMF_select;AMF_select];
                break
            end
        end
    end
    for j=2:size(raw,1)
        if(0<cell2mat(raw(j,16))&&cell2mat(raw(j,18))==0)
            UMF_select=[UMF_select;raw(j,:)];
        elseif(0<cell2mat(raw(j,16))+cell2mat(raw(j,18)))
            UMF_select1=[UMF_select1;raw(j,:)];
        end
    end
    
    
    cal_numberj=0; 
    A=ones(size(UMF_select,1),size(UMF_select1,1));
    B=[];
    C=zeros(2,size(UMF_select,1));
    for j=1:size(UMF_select,1)
        cal_numberk=0;
        for k=1:size(UMF_select1,1)
            A(j,k)=sum(abs(cell2mat(UMF_select(j,[6:11,20]))-cell2mat(UMF_select1(k,[6:11,20]))));
            if(A(j,k)<0.000001&&cell2mat(UMF_select(j,16))+cell2mat(UMF_select(j,18))==cell2mat(UMF_select1(k,16))+cell2mat(UMF_select1(k,18)))
                if(cell2mat(UMF_select(j,16))==1)
                    k1=1/Cl1_0;
                elseif(cell2mat(UMF_select(j,16))==2)
                    k1=1/Cl1_1;
                elseif(cell2mat(UMF_select(j,16))==3)
                    k1=1/Cl2_0;
                end
                if((cell2mat(UMF_select(j,1))/cell2mat(UMF_select1(k,1))<k1*(1+percentage_error)&&k1*(1-percentage_error)<cell2mat(UMF_select(j,1))/cell2mat(UMF_select1(k,1)))==1)
                    cal_numberk=cal_numberk+1;
                    C(cal_numberk,j)=k;
                end
            end
        end
        if(cal_numberk>0)
            cal_numberj=cal_numberj+1;
        end
        B(j)=cal_numberk;
    end
    UMF_result=cell(cal_numberj,1+2*max(B));
    
    cal_numberh=0;
    for j=1:length(B)
        if(B(j)>0)
            cal_numberh=cal_numberh+1;
            UMF_result{cal_numberh,1}=UMF_select(j,:);

            if(B(j)==1)
                UMF_result{cal_numberh,2}=[cell2mat(UMF_select1(C(1,j),[1:2,18])),cell2mat(UMF_select1(C(1,j),1))/cell2mat(UMF_select(j,1)),(cell2mat((UMF_select1(C(1,j),2)))-cell2mat(UMF_select(j,2))-deltam)/deltam*1000000];
                
            else
                UMF_result{cal_numberh,2}=[cell2mat(UMF_select1(C(1,j),[1:2,18])),cell2mat(UMF_select1(C(1,j),1))/cell2mat(UMF_select(j,1)),(cell2mat((UMF_select1(C(1,j),2)))-cell2mat(UMF_select(j,2))-deltam/deltam)*1000000];
                UMF_result{cal_numberh,4}=[cell2mat(UMF_select1(C(2,j),[1:2,18])),cell2mat(UMF_select1(C(1,j),1))/cell2mat(UMF_select(j,1)),(cell2mat((UMF_select1(C(1,j),2)))-cell2mat(UMF_select(j,2))-deltam/deltam)*1000000]; 
            end
        end
    end
    name_toprow={raw{1,:},'abundance','exp_mass','m_ppm','aband_Cl37/Cl','inten_error'};
    raw_judgediff=UMF_result;
    bottle1={};
    bottle2={};
    for j=1:size(raw_judgediff,1)
        cal_number=cal_number+1;
        bottle1=raw_judgediff{j,1};
        bottle2=raw_judgediff{j,2};
        diff_ppm(cal_number)=abs(abs(cell2mat(bottle1(2))-bottle2(2))-deltam)/cell2mat(bottle1(2))*1000000;
    end
end
diff_ppm=sort(diff_ppm);
ppm_allow1=max(20,diff_ppm(round(length(diff_ppm)*0.90)));
ppm_allow2=min(0,diff_ppm(round(length(diff_ppm)*0.1)));    
end

