function [ u0, count, far_field_value_u0, MIndices, I, J, b] = initial_approximation_Modificado1(h, R, MI); 

disp('Running: initial_approximation_Modificado1.m')

format short

%
%   Mesh parameters:
%

%h = 5e-2 % 
h= 0.5
%h = 5e-1
%N = 200
%R = N*h original
R=5
N = R/h;
x = -R:h:R  % x has 2N + 1 points

y = -R:h:R  % y has 2N + 1 points

NT = size(x,2) % NT = 2*N+1
n = 2;     % <--- space dimension is 2 (plane)
p = 4;     % <--- we are interested in the case: p > n
%format compact
%format long
gamma = (p - n)/(p - 1)   % <--- this parameter will not actually 
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

x1 = MI(1,1);   % <--- x-coordinate of 1st singularity
%x1= -9.95
y1 = MI(1,2);   % <--- y-coordinate of 1st singularity

i1 = find( x1-h/100 < x & x < x1 + h/100)  % <--- i-index of 1st singularity
j1 = find( y1-h/100 < y & y < y1 + h/100)   % <--- j-index of 1st singularity

x2 = MI(2,1);   % <--- x-coordinate of 2nd singularity
y2 = MI(2,2);  % <--- y-coordinate of 2nd singularity

i2 = find( x2-h/100 < x & x < x2 + h/100)   % <--- i-index of 2nd singularity
j2 = find( y2-h/100 < y & y < y2 + h/100)  % <--- j-index of 2nd singularity

x3 = MI(3,1);   % <--- x-coordinate of 3rd singularity
y3 = MI(3,2); % <--- y-coordinate of 3rd singularity

i3 = find( x3-h/100 < x & x < x3 + h/100)   % <--- i-index of 3rd singularity
j3 = find( y3-h/100 < y & y < y3 + h/100)  % <--- j-index of 3rd singularity

MIndices = [ i1 j1; i2 j2; i3 j3 ];
% -------------------------------------------------

b1 = MI(1,3);

b2 = MI(2,3);

b3 = MI(3,3);

%b = (b1 + b2 + b3)/3 -0.5% original
b = (b1 + b2 + b3)/3 - 1.5; 

% -------------------------------------------------

r12 = sqrt( (x1 - x2)^2 + (y1 - y2)^2 )
r13 = sqrt( (x1 - x3)^2 + (y1 - y3)^2 )
r23 = sqrt( (x2 - x3)^2 + (y2 - y3)^2 )

r0 = 1/2 * min([r12, r13, r23])

% -------------------------------------------------

u0 = b * ones(2*N+1,2*N+1);

for i = 1:2*N+1 % i é o x
    for j = 1:2*N+1 % j é o y
        % d1 = distãncia de cada ponto do grid a 1ª singularidade
        d1 = sqrt( (x(i) - x1)^2 + (y(j) - y1)^2 );
        
        d2 = sqrt( (x(i) - x2)^2 + (y(j) - y2)^2 );
        
        d3 = sqrt( (x(i) - x3)^2 + (y(j) - y3)^2 );
        
        if d1 < r0 % linha = i = x
           u0(i,j) = b1 - (b1 - b)*d1/r0;
        end
        
        if d2 < r0 % coluna = j = y
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

u0(I,J); % adicionado
M = [1 2 3; 4 5 6]; %Exemplo de matriz e exemplo de sum
M(:);
M(2,:); % 2ª linha
M(:,1); % 1ªcoluna
M(:,2); % 2ª coluna
size(M),size(M,1), sum(M,1), sum(M,2), sum(sum(M))
b;

count = 1;

min_u(count) = min( u0(:) );
max_u(count) = max( u0(:) );
min_u, max_u, size(min_u),size(max_u)
% -------------------------------------------------
%
%    Solution checkings:
%
variation_L1 (1) = NaN   % <--- L1 norm of u - u_previous
variation_sup(1) = NaN   % <--- supnorm of u - u_previous

%far_field_value(1) = far_field_value_u0

% -------------------------------------------------

%{ 
Example: 
Evaluate the function  x*exp(-x^2-y^2) 
over the range  -2 < x < 2,  -4 < y < 4,

figure(5)
[XX,YY] = meshgrid(-2:.2:2, -4:.4:4);
ZZ = XX .* exp(-XX.^2 - YY.^2);
surf(XX,YY,ZZ)
%}


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
u0;
% -------------------------------------------------
%{
v0 = u0;   % <--- this will be used by the code u_24h.m

t(1) = 0;  % <--- vector of time values for future solution checkings
t0 = 0;    % <--- this will be used by the code u_24h.m

dt_dump = 0.01;  % <--- time interval elapsed between solution checkings

% -------------------------------------------------
%}

 