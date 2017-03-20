%LAGG.M
%Y = lagg(y,L) where y is a time series of dimension Txv (with the first row 
%correspoinding to the earliest date) and L>=1 is an integer, produces a matrix 
%of order (T-L)x(v*(L+1)) where the first v columns are the contemporaneous time sereies
% the second v column are the original time series lagged 1 periods, and so on.
%The default value of L is 1.


function Y = lagg(y,L);

if nargin<2
L=1;
end

Y = [];
[T,v] = size(y);

for j=0:L
Y = [Y y(L+1-j:end-j,:)];
end

