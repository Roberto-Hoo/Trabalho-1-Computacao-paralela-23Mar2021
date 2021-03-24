function [T] = Tempo(m,n,p,R);
%
%
%
%

if (n == 2)
  Wn=2*pi;
elseif (n ==3)
    Wn =4*pi;
end %end if
%C = (m/Wn)

%C = (m/Wn) * (p/(p-1))

C = (m/Wn) * (p/(p-1))* ((p+n*(p-2))^(-n/p)) * ((p-2)/p)^(n*(p-1)/p);

C = C / beta(n*(p-1)/p , (2*p-3)/(p-2));

C = C ^ ((p*(p-2)) / ((p-1)*(p+n*(p-2))));

T = R * ( (p+n*(p-2)) ^ (-1/p) ) * ( ((p-2)/p) ^ ((p-1)/p) ) * ( C ^( 1/p -1));

T = T ^ ( p + n*(p-2) );

 