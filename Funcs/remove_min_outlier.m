function [range_processed] = remove_min_outlier(range, MIN_RANGE) 
% Function: remove outlier from LIDAR rang data.
%     range - range data.
%     MIN_RANGE - 0.03;
% Output:
%     range_processed - processed range data.
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), LI, Shuqing (jay_issaac@163.com), 20191202

if MIN_RANGE < 0.03
    MIN_RANGE = 0.03;
end
ref_size = 3;
w1 = 1.0/3; w2 = 1.0/3; w3 = 1.0/3; 
rang_rows = size(range,1); 
rang_cols = size(range,2);
% range_processed = zeros(size(range));
range_processed = zeros(rang_rows,rang_cols);
validRow = 1;
for j=1:rang_rows
    % search for valid points
    for i = 1:rang_cols-ref_size
        if range(j,i) > MIN_RANGE && range(j,i+1) > MIN_RANGE && range(j,i+2) > MIN_RANGE   % valid points
            break;
        end
    end
    validStart = i;
    if validStart >= rang_cols-ref_size % all data is invalid
        if validRow < j
            % copy last line
            range_processed(j,:) = range_processed(j-1,:);
        else
            % continue find valid lines
            validRow = validRow + 1;
            continue; 
        end
    end

    % fill outlier with left points
    for i = validStart:rang_cols
        if range(j,i) < MIN_RANGE    % outlier
            range_processed(j,i) = range_processed(j,i-1)*w1 + range_processed(j,i-2)*w2 + range_processed(j,i-3)*w3;
        else
            range_processed(j,i) = range(j,i);
       end
    end
    % fill outlier with right points and left points
    for i = rang_cols-ref_size :-1: validStart
        if range(j,i) < MIN_RANGE    % outlier
            range_processed(j,i) = (range_processed(j,i) + range_processed(j,i-1)*w1 + range_processed(j,i-2)*w2 + range_processed(j,i-3)*w3)/2;
        end
    end
    % fill outlier with right points at left edges
    for i = 1:validStart
        range_processed(j,i) = range_processed(j,i+1)*w1 + range_processed(j,i+2)*w2 + range_processed(j,i+3)*w3;
    end
end
for j=1:validRow-1
    range_processed(j,:) = range_processed(validRow,:);
end