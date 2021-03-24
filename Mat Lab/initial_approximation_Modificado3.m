clc % limpa tela
clear all

format compact
format long

%
%   Mesh parameters:
%

h = 0.1
R = 10
N = R/h

x = -R:h:R+h/100   % x has 2N + 1 points
y = -R:h:R+h/100;   % y has 2N + 1 points

n = 2;     % <--- space dimension is 2 (plane)
p = 4;     % <--- we are interested in the case: p > n

gamma = (p - n)/(p - 1);   % <--- this parameter will not actually 
                           %      be used in the program below
%
% -------------------------------------------------
%
%             INITIAL APPROXIMATION
%                     to
%             STEADY STATE SOLUTION
%
% -------------------------------------------------
%

% -------------------------------------------------

x1 = +1.0;   % <--- x-coordinate of 1st singularity
y1 = +0.5;   % <--- y-coordinate of 1st singularity

i1 = find( x1-h/100 < x & x < x1 + h/100);   % <--- i-index of 1st singularity
j1 = find( y1-h/100 < y & y < y1 + h/100);   % <--- j-index of 1st singularity

x2 = -1.5;   % <--- x-coordinate of 2nd singularity
y2 = +1.0;   % <--- y-coordinate of 2nd singularity

i2 = find( x2-h/100 < x & x < x2 + h/100);   % <--- i-index of 2nd singularity
j2 = find( y2-h/100 < y & y < y2 + h/100);   % <--- j-index of 2nd singularity

x3 = +0.5;   % <--- x-coordinate of 3rd singularity
y3 = -1.5;   % <--- y-coordinate of 3rd singularity

i3 = find( x3-h/100 < x & x < x3 + h/100);   % <--- i-index of 3rd singularity
j3 = find( y3-h/100 < y & y < y3 + h/100);   % <--- j-index of 3rd singularity

% -------------------------------------------------

b1 = 6;

b2 = 5;

b3 = 1;

b = (b1 + b2 + b3)/3 - 0.5;    

% -------------------------------------------------

r12 = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
r13 = sqrt( (x1 - x3)^2 + (y1 - y3)^2 );
r23 = sqrt( (x2 - x3)^2 + (y2 - y3)^2 );

r0 = 1/2 * min([r12, r13, r23])

% -------------------------------------------------

u0 = b * ones(2*N+1,2*N+1);

for i = 1:2*N+1
    for j = 1:2*N+1
        
        d1 = sqrt( (x(i) - x1)^2 + (y(j) - y1)^2 );
        
        d2 = sqrt( (x(i) - x2)^2 + (y(j) - y2)^2 );
        
        d3 = sqrt( (x(i) - x3)^2 + (y(j) - y3)^2 );
        
        if d1 < r0
           u0(i,j) = b1 - (b1 - b)*d1/r0;
        end
        
        if d2 < r0
           u0(i,j) = b2 - (b2 - b)*d2/r0;
        end  
        
        if d3 < r0
           u0(i,j) = b3 - (b3 - b)*d3/r0;
        end 
        
   end
    
end   
       
u0(i1,j1) = b1;
u0(i2,j2) = b2;
u0(i3,j3) = b3;

% -------------------------------------------------

I = find( 2.5 < abs(x) & abs(x) < 3.5 )
J = find( 2.5 < abs(y) & abs(y) < 3.5 )

far_field_value_u0 = sum( sum( u0(I,J) ) )/(length(I)*length(J))

b

count = 1;

min_u(count) = min( u0(:) );
max_u(count) = max( u0(:) );

% -------------------------------------------------
%
%    Solution checkings:
%
variation_L1 (1) = NaN;   % <--- L1 norm of u - u_previous
variation_sup(1) = NaN;   % <--- supnorm of u - u_previous

far_field_value(1) = far_field_value_u0;

% -------------------------------------------------

tic

[X,Y] = meshgrid(x,y);

figure(10)
mesh(X,Y,u0)
colormap(jet)
shading interp
xlabel('y')
ylabel('x')

figure(20)
shading interp
colormap gray
% surfl(X,Y,u0)
mesh(X,Y,u0)
xlabel('y')
ylabel('x')

toc

% -------------------------------------------------

v0 = u0;   % <--- this will be used by the code u_24h.m

t(1) = 0;  % <--- vector of time values for future solution checkings
t0 = 0;    % <--- this will be used by the code u_24h.m

dt_dump = 0.01;  % <--- time interval elapsed between solution checkings

% -------------------------------------------------
 m = pi * (r0^2) * abs(b1-b) * (1/3)
 n , p=4, R
[T] = Tempo(m,n,p,R)
%u0(107,86)
