function [rData_processed] = replace_outlier(rData, minRange, longRange) 
% Function: set outliers to a long range.
%     rData - range data (H*V).
%     minRange - minimal valid range (m).
%     longRange - outlier is set to this value (m).
% Output:
%     rData_processed - processed range data (H*V).
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200722
% 
[H,V] = size(rData);
rData_processed = rData;
for i = 1:H
    for j=2:V-1
        if rData_processed(i,j) < minRange
            rData_processed(i,j) = longRange;
        end
    end
end

