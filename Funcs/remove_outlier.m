function [rData_processed] = remove_outlier(rData, minRange) 
% Function: set outliers to neighboring values.
%     rData - range data (H*V).
%     minRange - minimal valid range (m).
% Output:
%     rData_processed - processed range data (H*V).
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200718
% 
[H,V] = size(rData);
rData_processed = rData;
for i = 1:H
    for j=2:V-1
        if rData_processed(i,j) < minRange
            % check previous
            if rData_processed(i,j-1) < minRange
                % check next
                if rData_processed(i,j+1) < minRange
                else
                    rData_processed(i,j) = rData_processed(i,j+1);
                end
            else
                % check next
                if rData_processed(i,j+1) < minRange
                    rData_processed(i,j) = rData_processed(i,j-1);
                else
                    % set to average
                    rData_processed(i,j) = ...
                        (rData_processed(i,j-1) + rData_processed(i,j-1))/2;
                end
            end
        end
    end
end

