function [range_processed] = replace_outlier(range, minRange, longRange) 
% Function: set outliers to a long range.
%     range - range data (H*V).
%     minRange - minimal valid range (m).
%     longRange - outlier is set to this value (m).
% Output:
%     range_processed - processed range data (H*V).
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200722
% 
[H,V] = size(range);
range_processed = range;
for i = 1:H
    for j=2:V-1
        if range_processed(i,j) < minRange
            range_processed(i,j) = longRange;
        end
    end
end

