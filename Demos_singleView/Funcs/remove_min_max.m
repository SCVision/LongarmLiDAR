function [range_processed] = remove_min_max(range, minRange, maxRange) 
% Function: set out-of-scope rang to 0.
%     range - range data.
%     minRange, maxRange - valid range scope
% Output:
%     range_processed - processed range data.
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200216
% Note: useless
if maxRange < 0
    maxRange = inf;
end
if minRange < 0
    minRange = 0;
end
inscope = (range < maxRange) & (range > minRange);
range_processed = range.*inscope;
