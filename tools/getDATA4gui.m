function [Y,X] = getDATA4gui(file,Ycol,Xcols)
% USAGE:
% - outputs:
%     Y: response vector;
%     X: predictors;
% - inputs:
%     file: name of the file where data is stored;
%     Ycol (OPTIONAL): specifies which column of the data array is to be
%     considered as the response. If no argument is given, the first column
%     is considered as the response.
%     Xcols (OPTIONAL): specifies which columns of the data array are to be
%     considered as predictors. If no argument is given, all columns aside
%     the response are considered as predictors.
%
% Notice this function is very close to the getDATA function. The only
% difference aims at reading formatted data regardless of the delimiter 
% used in the file. 
% =========================================================================

data = loadDATA(file);

if nargin < 3,
    Ycol = 1;
end
Yaux = data(:,Ycol);

if nargin < 4,
    data(:,Ycol) = [];
    X = data;
else
    X = data(:,Xcols);
end

Y = zeros(size(Yaux));
if isinZ(Yaux),
    minY = min(Yaux); maxY = max(Yaux);
    newidx = 1;
    for idx = minY:maxY,
        Y(Yaux==idx) = newidx;
        newidx = newidx+1;
    end
else
    Y = Yaux;
end

