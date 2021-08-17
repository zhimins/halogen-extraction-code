function [ordered_cell] = sortcell(raw_cell,row_start,sort_column,trend)
%Abstract: sort cell data by one of its columns with ascending or descending
%order
%   row_start: the part of the cell data preparing to sort, from which row to start cell
%   sort_column: use which column as order
%   trend: ascending or descending
unsorted_part=raw_cell(row_start:end,:);
blank=[];
switch trend
    case 'ascend'
        for i=1:size(unsorted_part,1)-1
            for j=(i+1):size(unsorted_part,1)
                if(unsorted_part{i,sort_column}>unsorted_part{j,sort_column})
                    blank=unsorted_part(j,:);
                    unsorted_part(j,:)=unsorted_part(i,:);
                    unsorted_part(i,:)=blank;
                else
                    continue
                end
            end
        end
    case 'decend'
        for i=1:size(unsorted_part,1)-1
            for j=(i+1):size(unsorted_part,1)
                if(unsorted_part{i,sort_column}<unsorted_part{j,sort_column})
                    blank=unsorted_part(j,:);
                    unsorted_part(j,:)=unsorted_part(i,:);
                    unsorted_part(i,:)=blank;
                else
                    continue
                end
            end
        end
end
ordered_cell=[raw_cell(1:(row_start-1),:);unsorted_part];
end

