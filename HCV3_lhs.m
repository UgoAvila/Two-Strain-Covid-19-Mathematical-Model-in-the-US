function  [X_scaled,X_normalized]=HCV3_lhs(n,min_ranges_p,max_ranges_p)
%	rangmin=[0.257 4.90195*10^8 0.02 10^4 0.01 0.009 10^2 10^(-4) 0.9 5*10^(-4) 17
%	2000 0.83178 10^(-8)]
%   rangmax=[1.028 9.9839*10^8 0.4 2*10^5 0.18 0.18 2*10^3 2*10^(-3) 18 4*10^(-3)
%   68 5000 16.6356 2*10^(-7)]
%   M=lhs(1000,rangmin,rangmax)
p=length(min_ranges_p);
[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(n,p);
OFFSET=ones(n,p);

for i=1:p
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
X_normalized = lhsdesign(n,p);

X_scaled=SLOPE.*X_normalized+OFFSET;