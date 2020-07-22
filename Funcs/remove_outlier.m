function [range_processed] = remove_outlier(range, minRange) 
% Function: set outliers to neighboring values.
%     range - range data (H*V).
%     minRange - minimal valid range (m).
% Output:
%     range_processed - processed range data (H*V).
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200718
% 
[H,V] = size(range);
range_processed = range;
for i = 1:H
    for j=2:V-1
        if range_processed(i,j) < minRange
            % check previous
            if range_processed(i,j-1) < minRange
                % check next
                if range_processed(i,j+1) < minRange
                else
                    range_processed(i,j) = range_processed(i,j+1);
                end
            else
                % check next
                if range_processed(i,j+1) < minRange
                    range_processed(i,j) = range_processed(i,j-1);
                else
                    % set to average
                    range_processed(i,j) = ...
                        (range_processed(i,j-1) + range_processed(i,j-1))/2;
                end
            end
        end
    end
end

