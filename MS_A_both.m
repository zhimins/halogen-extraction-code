function [list_MS] = MS_A_both(list_MS,i,foldername,X0_1_1,X1_1_1,X2_1_1,ppm_allow1_1,X0_1_2,X1_1_2,...
    X2_1_2,ppm_allow1_2,percentage_error,delta_m_1,delta_m_2,recal_ppm)

format LONG
X0_1=[X0_1_1,X0_1_2];
X1_1=[X1_1_1,X1_1_2];
X2_1=[X2_1_1,X2_1_2];
ppm_allow1=[ppm_allow1_1,ppm_allow1_2];
delta_m=[delta_m_1,delta_m_2];

%scan the list and do the isotopic double check

splitname = strtok(list_MS(i,1).name,'.');
[~,~,raw_MS]=xlsread([foldername,'\data\dataMS\align\blanksub\',list_MS(i,1).name]);
[~,~,raw_UMF]=xlsread([foldername,'\data\dataUMF\',[splitname(1:size(splitname,2)-2),'All_UMF.CSV']]);
raw_col=size(raw_UMF,2);
name_toprow={raw_UMF{1,1:raw_col},'abundance','exp_mass','mass ppm error','isotopic ratio','abundance error(%)'};
counterX=0;counterXplus2=0;
raw_selectXUMF=cell(1,raw_col);
raw_selectXplus2UMF=cell(1,raw_col);
%%divide data to X(Cl35/Br79£© and Xplus2(Cl37/Br81)
for j=2:size(raw_UMF,1)
    if((raw_UMF{j,16}>0&&raw_UMF{j,18}<1)||(raw_UMF{j,19}>0&&raw_UMF{j,20}<1))
        counterX=counterX+1;
        raw_selectXUMF(counterX,:)=raw_UMF(j,1:raw_col);
    elseif(raw_UMF{j,18}>0 || raw_UMF{j,20}>0)
        counterXplus2=counterXplus2+1;
        raw_selectXplus2UMF(counterXplus2,:)=raw_UMF(j,1:raw_col);
    end
end
%%add AMF data into consider
if (exist([foldername,'\data\dataAMF\',splitname(1:size(splitname,2)-2),'All_AMF.csv'],'file')==2)
    [~,~,raw_AMF]=xlsread([foldername,'\data\dataAMF\',[splitname(1:size(splitname,2)-2),'All_AMF.csv']]);
    if(size(raw_AMF,1)>1)
        for j=2:size(raw_AMF,1)
            if((raw_AMF{j,16}>0&&raw_AMF{j,18}<1)||(raw_AMF{j,19}>0&&raw_AMF{j,20}<1))
                counterX=counterX+1;
                raw_selectXUMF(counterX,:)=raw_AMF(j,1:raw_col);
            elseif(raw_AMF{j,18}>0 || raw_AMF{j,20}>0)
                counterXplus2=counterXplus2+1;
                raw_selectXplus2UMF(counterXplus2,:)=raw_AMF(j,1:raw_col);
            end
        end
    end
end
% X-DBPs determined by isotope double check 
A_result_direct=[];
count_UMF_direct=0;
for j=size(raw_selectXUMF,1):-1:1
    for k=size(raw_selectXplus2UMF,1):-1:1
        if(strcmp(raw_selectXUMF{j,3},raw_selectXplus2UMF{k,3}) &&...
                ((raw_selectXplus2UMF{k,2}-raw_selectXUMF{j,2}-1.997)/raw_selectXUMF{j,2}*1000000<ppm_allow1(1) ||...
                (raw_selectXplus2UMF{k,2}-raw_selectXUMF{j,2}-1.998)/raw_selectXUMF{j,2}*1000000<ppm_allow1(2)) &&...
                ( ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X0_1(2)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X0_1(2)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,19))==1)) || ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X1_1(2)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X1_1(2)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,19))==2)) || ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X2_1(2)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X2_1(2)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,19))==3)) || ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X0_1(1)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X0_1(1)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,16))==1)) || ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X1_1(1)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X1_1(1)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,16))==2)) || ((raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}>X2_1(1)*0.7 && raw_selectXplus2UMF{k,1}/raw_selectXUMF{j,1}<X2_1(1)*1.3 &&...
                        cell2mat(raw_selectXUMF(j,16))==3))))
            count_UMF_direct=count_UMF_direct+1;
            A_result_direct{count_UMF_direct,1}=raw_selectXUMF(j,:);
            A_result_direct{count_UMF_direct,2}=[raw_selectXplus2UMF{k,1},raw_selectXplus2UMF{k,2},0,0,0];
            raw_selectXUMF(j,:)=[];
            raw_selectXplus2UMF(k,:)=[];
            break
        end
    end
end

%%search the corresponding isotopic mass in MS
cal_number=0;
UMF2MS={};
UMF2MS_mi={};
for j=1:size(raw_selectXUMF,1)
    for k=2:size(raw_MS,1)
        if(cell2mat(raw_MS(k,2))==cell2mat(raw_selectXUMF(j,1)) && ...
                abs(raw_MS{k,1}-raw_selectXUMF{j,2})/raw_selectXUMF{j,2}*1000000<recal_ppm)
            cal_number=cal_number+1;
            UMF2MS(cal_number,1:size(raw_selectXUMF,2))=raw_selectXUMF(j,:);
            UMF2MS_mi(cal_number,1:size(raw_MS,2))=raw_MS(k,:);
        end
    end
end

UMF2MS=sortcell(UMF2MS,1,2,'ascend');
UMF2MS_mi=sortcell(UMF2MS_mi,1,1,'ascend');

%%X-DBPs determined by isotope double check (Cl35/Br79)(mass scaning)
cal_number=0;cal_number_Cl=0;cal_number_Br=0;cal_number_ClBr=0;
A_select=[];A_select_Cl=[];A_select_Br=[];A_select_ClBr=[];
for j=1:size(UMF2MS,1)
    Limit=find(cell2mat(raw_MS(2:end,1))>(UMF2MS{j,2}-0.5));%only one Cl37 existed
    Limit2=find(cell2mat(raw_MS(2:end,1))<(UMF2MS{j,2}+2.5));
    Limit=sort(Limit);Limit2=sort(Limit2);
    if(isempty(Limit)==0)
        if(Limit(1)<=Limit2(end))
            for k=max(2,Limit(1)):min(Limit2(end),size(raw_MS,1))
                if(abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(1))/cell2mat(UMF2MS_mi(j,1))*1000000<=ppm_allow1(1) &&...
                        ( ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X0_1(1)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X0_1(1)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,16))==1 && cell2mat(UMF2MS(j,19))==0)) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X1_1(1)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X1_1(1)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,16))==2 && cell2mat(UMF2MS(j,19))==0)) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X2_1(1)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X2_1(1)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,16))==3 && cell2mat(UMF2MS(j,19))==0)) ) )
                    cal_number_Cl=cal_number_Cl+1;
                    A_select_Cl(1,cal_number_Cl)=j;
                    A_select_Cl(2,cal_number_Cl)=k;
                    A_select_Cl(3,cal_number_Cl)=abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(1))/cell2mat(UMF2MS_mi(j,1))*1000000;
                end
                if(abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(2))/cell2mat(UMF2MS_mi(j,1))*1000000<=ppm_allow1(2) &&...
                        ( ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X0_1(2)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X0_1(2)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,19))==1 && cell2mat(UMF2MS(j,16))==0)) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X1_1(2)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X1_1(2)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,19))==2 && cell2mat(UMF2MS(j,16))==0)) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=X2_1(2)*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=X2_1(2)*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,19))==3 && cell2mat(UMF2MS(j,16))==0)) ) )
                    cal_number_Br=cal_number_Br+1;
                    A_select_Br(1,cal_number_Br)=j;
                    A_select_Br(2,cal_number_Br)=k;
                    A_select_Br(3,cal_number_Br)=abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(2))/cell2mat(UMF2MS_mi(j,1))*1000000;
                end
                if(abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(1))/cell2mat(UMF2MS_mi(j,1))*1000000<=ppm_allow1(2) &&...
                        ( ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=1.2924*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=1.2924*(1+percentage_error)) &&...
                        ( cell2mat(UMF2MS(j,19))==1 && cell2mat(UMF2MS(j,16))==1 )) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=2.2652*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=2.2652*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,19))==2 && cell2mat(UMF2MS(j,16))==1 )) || ((cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))>=1.612*(1-percentage_error) && cell2mat(raw_MS(k,2))/cell2mat(UMF2MS_mi(j,2))<=1.612*(1+percentage_error) &&...
                        cell2mat(UMF2MS(j,19))==1 && cell2mat(UMF2MS(j,16))==2 )) ) )
                    cal_number_ClBr=cal_number_ClBr+1;
                    A_select_ClBr(1,cal_number_ClBr)=j;
                    A_select_ClBr(2,cal_number_ClBr)=k;
                    A_select_ClBr(3,cal_number_ClBr)=abs(cell2mat(raw_MS(k,1))-cell2mat(UMF2MS_mi(j,1))-delta_m(1))/cell2mat(UMF2MS_mi(j,1))*1000000;
                end
            end
        end
    end
end

A_select=[A_select_Cl,A_select_Br,A_select_ClBr];
A_result=cell(size(A_select,2),2);
for j=1:size(A_select,2)
    A_result{j,1}=UMF2MS(A_select(1,j),:);
    A_result{j,2}=[raw_MS(A_select(2,j),2),raw_MS(A_select(2,j),1),A_select(3,j)];
end
im=[];
%judge means X,judge1 means Xplus2
if(isempty(A_result)==0)
    for j=1:size(A_result,1)
        inten_judge=A_result{j,1};
        inten_judge1=A_result{j,2};
        Inten_judge=inten_judge{1};
        Inten_judge1=inten_judge1{1};
        im(j,1)=Inten_judge1(1)/Inten_judge(1);
    end
    
    inten_judge_=[];
    cal_numberjudge=0;
    intencal_number=0;
    intencal=[];
    
    for j=1:cal_number_Cl
        Inten_judge=A_result{j,1};
        if(Inten_judge{16}==1 && Inten_judge{18}==0 && Inten_judge{19}==0 && Inten_judge{20}==0)
            if(abs(im(j)-X0_1(1))/X0_1(1)<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X0_1(1);
                intencal(intencal_number,1)=abs(im(j)-X0_1(1))/X0_1(1)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{16}==2 && Inten_judge{18}==0 && Inten_judge{19}==0 && Inten_judge{20}==0)
            if(abs(im(j)-X1_1(1))/(X1_1(1))<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X1_1(1);
                intencal(intencal_number,1)=abs(im(j)-X1_1(1))/X1_1(1)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{16}==3 && Inten_judge{18}==0 && Inten_judge{19}==0 && Inten_judge{20}==0)
            if(abs(im(j)-X2_1(1))/(X2_1(1))<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X2_1(1);
                intencal(intencal_number,1)=abs(im(j)-X2_1(1))/X2_1(1)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        else
            cal_numberjudge=cal_numberjudge+1;
            inten_judge_(cal_numberjudge)=j;
            continue
        end
    end
    for j=(cal_number_Cl+1):(cal_number_Cl+cal_number_Br)%size(im,1)
        Inten_judge=A_result{j,1};
        if(Inten_judge{19}==1 && Inten_judge{20}==0 && Inten_judge{16}==0 && Inten_judge{18}==0)
            if(abs(im(j)-X0_1(2))/X0_1(2)<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X0_1(2);
                intencal(intencal_number,1)=abs(im(j)-X0_1(2))/X0_1(2)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{19}==2 && Inten_judge{20}==0 && Inten_judge{16}==0 && Inten_judge{18}==0)
            if(abs(im(j)-X1_1(2))/(X1_1(2))<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X1_1(2);
                intencal(intencal_number,1)=abs(im(j)-X1_1(2))/X1_1(2)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{19}==3 && Inten_judge{20}==0 && Inten_judge{16}==0 && Inten_judge{18}==0)
            if(abs(im(j)-X2_1(2))/(X2_1(2))<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=X2_1(2);
                intencal(intencal_number,1)=abs(im(j)-X2_1(2))/X2_1(2)*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        else
            cal_numberjudge=cal_numberjudge+1;
            inten_judge_(cal_numberjudge)=j;
            continue
        end
    end
    
    for j=(cal_number_Cl+cal_number_Br+1):size(im,1)
        Inten_judge=A_result{j,1};
        if(Inten_judge{19}==1 && Inten_judge{20}==0 && Inten_judge{16}==1 && Inten_judge{18}==0)
            if(abs(im(j)-1.2924)/1.2924<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=1.2924;
                intencal(intencal_number,1)=abs(im(j)-1.2924)/1.2924*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{19}==2 && Inten_judge{20}==0 && Inten_judge{16}==1 && Inten_judge{18}==0)
            if(abs(im(j)-2.2652)/2.2652<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=2.2652;
                intencal(intencal_number,1)=abs(im(j)-2.2652)/2.2652*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        elseif(Inten_judge{19}==1 && Inten_judge{20}==0 && Inten_judge{16}==2 && Inten_judge{18}==0)
            if(abs(im(j)-1.612)/1.612<=percentage_error)
                intencal_number=intencal_number+1;
                Abundance_ratio(intencal_number,1)=1.612;
                intencal(intencal_number,1)=abs(im(j)-1.612)/1.612*100;
                continue
            else
                cal_numberjudge=cal_numberjudge+1;
                inten_judge_(cal_numberjudge)=j;
                continue
            end
        else
            cal_numberjudge=cal_numberjudge+1;
            inten_judge_(cal_numberjudge)=j;
            continue
        end
    end
    
    for j=size(inten_judge_,2):-1:1
        for k=size(A_result,1):-1:1
            if(k==inten_judge_(j))
                A_result(k,:)=[];
            end
        end
    end
    
    
    for j=1:length(intencal)
        A_result{j,2}=[A_result{j,2},Abundance_ratio(j),intencal(j)];
    end
    xlswrite([foldername,'\output\MS_A\',splitname],name_toprow,1,'A1');
    A_result=[A_result_direct;A_result];
    
    Inten_judge={};
    Inten_judge1={};
    cal_number=0;
    for j=1:size(A_result,1)
        if(isempty(A_result{j,1})==0)
            if(cal_number>=1)
                if(strcmp(A_result{j,1}{3},Inten_judge{3})==0)
                    cal_number=cal_number+1;
                    Inten_judge=A_result{j,1};
                    Inten_judge1=A_result{j,2};
                    xlswrite([foldername,'\output\MS_A\',splitname],Inten_judge,1,['A',num2str(cal_number+1)]);
                    xlswrite([foldername,'\output\MS_A\',splitname],Inten_judge1,1,['AX',num2str(cal_number+1)]);
                end
            else
                cal_number=cal_number+1;
                Inten_judge=A_result{j,1};
                Inten_judge1=A_result{j,2};
                xlswrite([foldername,'\output\MS_A\',splitname],Inten_judge,1,['A',num2str(cal_number+1)]);
                xlswrite([foldername,'\output\MS_A\',splitname],Inten_judge1,1,['AX',num2str(cal_number+1)]);
            end
        end
    end
end
    
end

