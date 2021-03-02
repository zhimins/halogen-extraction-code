function [list_UMF]=UMF_unsureCl(list_MS,foldername)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明


format LONG
%list_MS=dir([foldername,'\output\MS_A\','*.xls']);
list_UMF=dir([foldername,'\data\dataUMF\','*.csv']);
for i=1:size(list_UMF,1)
    splitname = strtok(list_UMF(i,1).name,'.');
    [~,~,raw]=xlsread([foldername,'\data\dataUMF\',list_UMF(i,1).name]);
    [~,~,raw_right]=xlsread([foldername,'\output\MS_A\',[splitname(1:length(splitname)-7),'MS.xls']]);
    %raw_select=xlsread([foldername,'\select\',[splitname,'.xls']]);
    raw_select=cell(1,size(raw,2));
    raw_select(1,:)=raw(1,:);
    count=1;
    for j=2:size(raw,1)
        if(raw{j,16}+raw{j,18}>0)
            count=count+1;
            raw_select(count,:)=raw(j,:);
        end
    end
    for j=size(raw_right,1):-1:2
        for k=size(raw_select,1):-1:2
            if(raw_select{k,1}==raw_right{j,1})
                raw_select(k,:)=[];
            end
        end
    end
    for j=size(raw_select,1):-1:2
        for k=size(raw,1):-1:2
            if(raw{k,1}==raw_select{j,1})
                raw(k,:)=[];
            end
        end
    end
    %writecell(raw,[foldername,'\output\UMF\',splitname,'.csv']);% export the 
    %writecell(raw_select,[foldername,'\output\unsureCl\',[splitname,'_unsureCl.csv']]);
    cell2csv([foldername,'\output\UMF\',splitname,'.csv'], raw);
    cell2csv([foldername,'\output\unsureCl\',[splitname,'_unsureCl.csv']], raw_select);
end
end

