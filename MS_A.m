function [list_MS] = MS_A(foldername,Cl1_0,Cl1_1,Cl2_0,ppm_allow1,ppm_allow2,percentage_error,deltam)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
format LONG
list_MS=dir([foldername,'\data\dataMS\','*.csv']);%'*.CSV''*.csv'
h=waitbar(0,'MS_A is calculating')
for i=1:size(list_MS,1)
    
    splitname = strtok(list_MS(i,1).name,'.');
    [~,~,raw_MS]=xlsread([foldername,'\data\dataMS\',list_MS(i,1).name]);
    [~,~,raw_UMF]=xlsread([foldername,'\data\dataUMF\',[splitname(1:size(splitname,2)-2),'All_UMF.csv']]);
    raw_col=size(raw_UMF,2);
    name_toprow={raw_UMF{1,1:raw_col},'abundance','exp_mass','m_ppm','aband_Cl37/Cl','inten_error'};
    conter35=0;
    raw_select35UMF=cell(1,raw_col);%size(raw_UMF,2));
    for j=2:size(raw_UMF,1)
        if(raw_UMF{j,16}>0&&raw_UMF{j,18}<1)
            conter35=conter35+1;
            raw_select35UMF(conter35,:)=raw_UMF(j,1:raw_col);
        end
    end
    if (exist([foldername,'\data\dataAMF\',splitname(1:size(splitname,2)-2),'All_AMF.csv'],'file')==2)
        [~,~,raw_AMF]=xlsread([foldername,'\data\dataAMF\',[splitname(1:size(splitname,2)-2),'All_AMF.csv']]);
        if(size(raw_AMF,1)>1)
            for j=2:size(raw_AMF,1)
                if(raw_AMF{j,16}>0&&raw_AMF{j,18}<1)
                    conter35=conter35+1;
                    raw_select35UMF(conter35,:)=raw_AMF(j,1:raw_col);
                end
            end
        end
	    %mkdir('路径');
    end
    cal_number=0;
    UMF2MS={};
    UMF2MS_mi={};
    for j=1:size(raw_select35UMF,1)
        for k=2:size(raw_MS,1)
            if(cell2mat(raw_MS(k,2))==cell2mat(raw_select35UMF(j,1)))
                cal_number=cal_number+1;
                UMF2MS(cal_number,1:size(raw_select35UMF,2))=raw_select35UMF(j,:);
                UMF2MS_mi(cal_number,1:size(raw_MS,2))=raw_MS(k,:);
                break
            end
        end
    end
    cal_number=0;
    A_select=[];
    for j=1:size(UMF2MS,1)
        for k=2:size(raw_MS,1)
            if(abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-deltam)/cell2mat(UMF2MS_mi(j,1))*1000000<=ppm_allow1&&ppm_allow2<=abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-deltam)/cell2mat(UMF2MS_mi(j,1))*1000000)
                cal_number=cal_number+1;
                A_select(1,cal_number)=j;
                A_select(2,cal_number)=k;
                A_select(3,cal_number)=abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-deltam)/cell2mat(UMF2MS_mi(j,1))*1000000;
            end
        end
    end
    A_result=cell(size(A_select,2),2);
    for j=1:size(A_select,2)
        A_result{j,1}=UMF2MS(A_select(1,j),:);
        A_result{j,2}=[raw_MS(A_select(2,j),2),raw_MS(A_select(2,j),1),A_select(3,j)];
    end
    im=[];
    if(isempty(A_result)==0)
        for j=1:size(A_result,1)
            inten_judge=A_result{j,1};
            inten_judge1=A_result{j,2};
            Inten_judge=inten_judge{1};
            Inten_judge1=inten_judge1{1};
            im(j,1)=(Inten_judge1(1)/Inten_judge(1));
        end
        inten_judge_=[];
        cal_numberjudge=0;
        intencal_number=0;
        intencal=[];
        for j=1:size(im,1)
            Inten_judge=A_result{j,1};
            if(Inten_judge{18}+Inten_judge{16}==1)
                if(abs(im(j)-Cl1_0)/Cl1_0<percentage_error)
                    intencal_number=intencal_number+1;
                    intencal(intencal_number,1)=abs(im(j)-Cl1_0)/Cl1_0*100;
                    continue
                else
                    cal_numberjudge=cal_numberjudge+1;
                    inten_judge_(cal_numberjudge)=j;
                    continue
                end
            elseif(Inten_judge{18}+Inten_judge{16}==2)
                if(abs(im(j)-Cl1_1)/(Cl1_1)<percentage_error)
                    intencal_number=intencal_number+1;
                    intencal(intencal_number,1)=abs(im(j)-Cl1_1)/Cl1_1*100;
                    continue
                else
                    cal_numberjudge=cal_numberjudge+1;
                    inten_judge_(cal_numberjudge)=j;
                    continue
                end
            elseif(Inten_judge{18}+Inten_judge{16}==3)
                if(abs(im(j)-Cl2_0)/(Cl2_0)<percentage_error)
                    intencal_number=intencal_number+1;
                    intencal(intencal_number,1)=abs(im(j)-Cl2_0)/Cl2_0*100;
                    continue
                else
                    cal_numberjudge=cal_numberjudge+1;
                    inten_judge_(cal_numberjudge)=j;
                    continue
                end
            end
        end
        for j=size(inten_judge_,2):-1:1
            for k=size(A_result,1):-1:1
                if(k==inten_judge_(j))
                    A_result{k,1}={};
                    A_result{k,2}={};
                end
            end
        end
        xlswrite([foldername,'\output\MS_A\',splitname],name_toprow,1,'A1');
        Inten_judge={};
        Inten_judge1={};
        cal_number=0;
        for j=1:size(A_result,1)
            if(isempty(A_result{j,1})==0)
                if(cal_number>=1)
                    if(length(A_result{j,1}{3})~=length(Inten_judge{3}))
                        cal_number=cal_number+1;
                       Inten_judge=A_result{j,1};
                        Inten_judge1=A_result{j,2};
                        xlswrite([foldername,'\output\MS_A\',splitname],[Inten_judge,Inten_judge1,cell2mat(Inten_judge1(1))/cell2mat(Inten_judge(1)),intencal(cal_number)],1,['A',num2str(cal_number+1)]);
                    elseif(A_result{j,1}{3}~=Inten_judge{3})
                        cal_number=cal_number+1;
                       Inten_judge=A_result{j,1};
                        Inten_judge1=A_result{j,2};
                        xlswrite([foldername,'\output\MS_A\',splitname],[Inten_judge,Inten_judge1,cell2mat(Inten_judge1(1))/cell2mat(Inten_judge(1)),intencal(cal_number)],1,['A',num2str(cal_number+1)]);
                    end
                else
                     cal_number=cal_number+1;
                     Inten_judge=A_result{j,1};
                     Inten_judge1=A_result{j,2};
                     xlswrite([foldername,'\output\MS_A\',splitname],[Inten_judge,Inten_judge1,cell2mat(Inten_judge1(1))/cell2mat(Inten_judge(1)),intencal(cal_number)],1,['A',num2str(cal_number+1)]);
                end
            end
        end
    end
    waitbar(i/size(list_MS,1))
end
close(h)
end

