function [data4,data5]=merge_data(data1,data2) 
% Function: Merge pre-scanning data and adaptive scanning data.
% Input:
%     data1,data2 - Pre-scanning data and adaptive scanning data. (data=[angleH,range])
% Output:
%     data4 - Merged data. (data=[angleH,range])
%     data5 - The same angleH data in pre-scanning data and adaptive scanning data (X,Y,Z). (data=[angleH,range]) 
% Demo:
% [range1, angleV1, angleH1, timestamp1] = read_scandata("batchScanned20191126210113.txt");
% [range2, angleV2, angleH2, timestamp2] = read_scandata("batchScanned20191126212758.txt");
% data1=[angleH1,range1];
% data2=[angleH2,range2];
% [data4,data5]=merge_data(data1,data2); 
% 
% Writen by LI, Shuqing (jay_issaac@163.com), 20191205
%

    temp=[data1;data2];
    %% merge data and reorder
    [data3,I]=sortrows(temp,1,'descend');
    %% get mean valve of the same data
    j=1;
    k=1;
    l=1;
    [n3,m3] = size(data3);
    for i=1:n3
        if abs(data3(l+1,1)-data3(l,1))<0.000001
            data4(j,1:m3)=(data3(l,:)+data3(l+1,:))/2;
            data5(k,1:m3)=data3(l,:);
            data5(k+1,1:m3)=data3(l+1,:);
            k=k+2;
            l=l+2;
            if l>n3
                break;
            end
            j=j+1;
        else
            data4(j,1:m3)=data3(l,:);
            j=j+1;
            l=l+1;
            if l>n3
                break;
            end
        end
    end
end