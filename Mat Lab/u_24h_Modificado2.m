
%
% **************************************************
%
%    UPDATE LINES   14-15 183-260
%
% **************************************************
%
%  Code INITIAL_APPROXIMATION must have already been run !!! <---IMPORTANT
%
% **************************************************
%

% ************************************************** <--|
%                                                       |
% load data_tF47p0      % <--- needs updating           | please ignore
% u0 = u_47p0;          % <--- needs updating           | these steps
%                                                       | (for now)
% ************************************************** <--|

%
% Computing steady-state limit of solution u(.,t):
%

cfl = 0.01;   % <--- Courant-Friedrichs-Lewy number

dt = cfl*h*h;  % <--- time step length

tF = t0 + 1; 

% --------------------------------------------------------------------

no_time_steps = round((tF-t0)/dt);

iter_max = round((tF-t0)/dt_dump);

iterations_per_cycle = round(dt_dump/dt);

format compact
format long

disp(' ')
disp('*******************************************************************')
disp('Computing ...')
disp('*******************************************************************')
disp(' ')
disp('---------------------------------------------------------------')
disp('    variation in L1     variation (sup)     far field value    ')
disp('---------------------------------------------------------------')
disp([variation_L1(1:count), variation_sup(1:count), far_field_value(1:count)]) 
disp('---------------------------------------------------------------')

u = v0;
u_previous = v0;

q = (p - 2)/2;
 

for iter = 1:iter_max
  
  tic
  
  count = count + 1;
  
  t(count) = t(count - 1) + iterations_per_cycle*dt;

  v = zeros(2*N+1,2*N+1);

  F = zeros(2*N+1,2*N+1);

  G = zeros(2*N+1,2*N+1);
 
  
  for k = 1:iterations_per_cycle
 
  
    % ****************************************************************** %
    %                computation of F values on the grid:                %
    % ****************************************************************** %
  
    % ******  F(i,:) for 2 <= i <= 2*N + 1:
  
      F(2:2*N+1,2:2*N) = ( (u(2:2*N+1,2:2*N) - u(1:2*N,2:2*N)).^2 + ...
      ( (u(2:2*N+1,3:2*N+1)+u(1:2*N,3:2*N+1))-(u(2:2*N+1,1:2*N-1)+u(1:2*N,1:2*N-1)) ).^2 /16 ).^q / h^(p-2);
                   %  for 2 <= j <= 2*N
                 
      F(2:2*N+1,1) = ( (u(2:2*N+1,1) - u(1:2*N,1)).^2 + ...
      ( (u(2:2*N+1,2)+u(1:2*N,2))-(u(2:2*N+1,1)+u(1:2*N,1)) ).^2 /16 ).^q / h^(p-2);
                   %  for j = 1
                 
      F(2:2*N+1,2*N+1) = ( (u(2:2*N+1,2*N+1) - u(1:2*N,2*N+1)).^2 + ...
      ( (u(2:2*N+1,2*N+1)+u(1:2*N,2*N+1))-(u(2:2*N+1,2*N)+u(1:2*N,2*N)) ).^2 /16 ).^q / h^(p-2);
                   %  for j = 2*N + 1
                 
    % ******  F(i,:) for i = 1:
  
      F(1,2:2*N) = ( (u(1,3:2*N+1) - u(1,1:2*N-1)).^2 /4 ).^q / h^(p-2);
                   %  for 2 <= j <= 2*N
                 
      F(1,1) = ( (u(1,2) - u(1,1))^2 / 4 )^q / h^(p-2);
                   %  for j = 1
                 
      F(1,2*N+1) = ( (u(1,2*N+1) - u(1,2*N))^2 / 4 )^q / h^(p-2);
                   %  for j = 2*N + 1
  
                 
    % ****************************************************************** %
    %                computation of G values on the grid:                %
    % ****************************************************************** %
  
    % ******  G(:,j) for 2 <= j <= 2*N + 1: 
    
      G(2:2*N,2:2*N+1) = ( (u(2:2*N,2:2*N+1) - u(2:2*N,1:2*N)).^2 + ...
      ( (u(3:2*N+1,2:2*N+1)+u(3:2*N+1,1:2*N))-(u(1:2*N-1,2:2*N+1)+u(1:2*N-1,1:2*N)) ).^2 /16 ).^q / h^(p-2);     
                   %  for 2 <= i <= 2*N
                 
      G(1,2:2*N+1) = ( (u(1,2:2*N+1) - u(1,1:2*N)).^2 + ...
      ( (u(2,2:2*N+1)+u(2,1:2*N))-(u(1,2:2*N+1)+u(1,1:2*N)) ).^2 /16 ).^q / h^(p-2); 
                   %  for i = 1
    
      G(2*N+1,2:2*N+1) = ( (u(2*N+1,2:2*N+1) - u(2*N+1,1:2*N)).^2 + ...
      ( (u(2*N+1,2:2*N+1)+u(2*N+1,1:2*N))-(u(2*N,2:2*N+1)+u(2*N,1:2*N)) ).^2 /16 ).^q / h^(p-2);     
                   %  for i = 2*N + 1
    
    % ******  G(:,j) for j = 1:
                 
      G(2:2*N,1) = ( (u(3:2*N+1,1) - u(1:2*N-1,1)).^2 /4 ).^q ./ h^(p-2);            
                   %  for 2 <= i <= 2*N
  
      G(1,1) = ( (u(2,1) - u(1,1))^2 /4 )^q / h^(p-2);    
                   %  for i = 1
  
      G(2*N+1,1) = ( (u(2*N+1,1) - u(2*N,1)).^2 /4 )^q / h^(p-2);            
                   %  for 2 <= i <= 2*N
  
                 
    % ****************************************************************** %
    %      computation of v = [ new u values at the new time level ]:    %
    % ****************************************************************** %
  
    % ******  v(i,:) for 2 <= i <= 2*N: 
    
      v(2:2*N,2:2*N) = u(2:2*N,2:2*N) + ...
          cfl*( F(3:2*N+1,2:2*N).*(u(3:2*N+1,2:2*N)-u(2:2*N,2:2*N)) - ...    
                F(2:2*N,2:2*N).*(u(2:2*N,2:2*N)-u(1:2*N-1,2:2*N)) ) + ...
          cfl*( G(2:2*N,3:2*N+1).*(u(2:2*N,3:2*N+1)-u(2:2*N,2:2*N)) - ...
                G(2:2*N,2:2*N).*(u(2:2*N,2:2*N)-u(2:2*N,1:2*N-1)) ); 
                    %  for 2 <= j <= 2*N
                
      v(2:2*N,1) = u(2:2*N,1) + ...
          cfl*( F(3:2*N+1,1).*(u(3:2*N+1,1)-u(2:2*N,1)) - ...    
                F(2:2*N,1).*(u(2:2*N,1)-u(1:2*N-1,1)) ) + ...                 
          cfl*( G(2:2*N,2).*(u(2:2*N,2)-u(2:2*N,1)) );
                    %  for j = 1
                
      v(2:2*N,2*N+1) = u(2:2*N,2*N+1) + ...
          cfl*( F(3:2*N+1,2*N+1).*(u(3:2*N+1,2*N+1)-u(2:2*N,2*N+1)) - ...    
                F(2:2*N,2*N+1).*(u(2:2*N,2*N+1)-u(1:2*N-1,2*N+1)) ) - ...                 
          cfl*( G(2:2*N,2*N+1).*(u(2:2*N,2*N+1)-u(2:2*N,2*N)) );
                    %  for j = 2*N + 1
   
    % ******  v(i,:) for i = 1: 
    
      v(1,2:2*N) = u(1,2:2*N) + ...
          cfl*( F(2,2:2*N).*(u(2,2:2*N)-u(1,2:2*N)) ) + ...    
          cfl*( G(1,3:2*N+1).*(u(1,3:2*N+1)-u(1,2:2*N)) - ...
                G(1,2:2*N).*(u(1,2:2*N)-u(1,1:2*N-1)) ); 
                    %  for 2 <= j <= 2*N            
                
      v(1,1) = u(1,1) + ...
          cfl*( F(2,1)*(u(2,1)-u(1,1)) ) + ...    
          cfl*( G(1,2)*(u(1,2)-u(1,1)) );
                    %  for j = 1       
    
      v(1,2*N+1) = u(1,2*N+1) + ...
          cfl*( F(2,2*N+1)*(u(2,2*N+1)-u(1,2*N+1)) ) - ...    
          cfl*( G(1,2*N+1)*(u(1,2*N+1)-u(1,2*N)) );
                    %  for j = 2*N + 1
                               
    % ******  v(i,:) for i = 2*N + 1: 
    
      v(2*N+1,2:2*N) = u(2*N+1,2:2*N) + ...
          cfl*( - F(2*N+1,2:2*N).*(u(2*N+1,2:2*N)-u(2*N,2:2*N)) ) + ...    
          cfl*( G(2*N+1,3:2*N+1).*(u(2*N+1,3:2*N+1)-u(2*N+1,2:2*N)) - ...
                G(2*N+1,2:2*N).*(u(2*N+1,2:2*N)-u(2*N+1,1:2*N-1)) ); 
                    %  for 2 <= j <= 2*N            
                
      v(2*N+1,1) = u(2*N+1,1) + ...
          cfl*( - F(2*N+1,1)*(u(2*N+1,1)-u(2*N,1)) ) + ...    
          cfl*( G(2*N+1,2)*(u(2*N+1,2)-u(2*N+1,1)) );
                    %  for j = 1       
    
      v(2*N+1,2*N+1) = u(2*N+1,2*N+1) + ...
          cfl*( - F(2*N+1,2*N+1)*(u(2*N+1,2*N+1)-u(2*N,2*N+1)) ) - ...    
          cfl*( G(2*N+1,2*N+1)*(u(2*N+1,2*N+1)-u(2*N+1,2*N)) );
                    %  for j = 2*N + 1               
 
                    
   % ****************************************************************** %
   %            computation of new u values completed!                  %
   % ****************************************************************** %
                              
      u = v;    % <--- updating u (at the new time level)
    
      % correcting u values at the singular points:
    
      u(i1,j1) = b1;
    
      u(i2,j2) = b2;
    
      u(i3,j3) = b3;
  
   end  

   time_per_cycle(iter) = toc;
   
   % computing far field values of current solution:
   
   far_field_value(count) = sum( sum( u(I,J) ) )/(length(I)*length(J));
   
   % saving minimum and maximum solution values at current time:
   
   min_u(count) = min( u(:) );
   max_u(count) = max( u(:) );
 
   % 
   mass_u(count) = sum( u(:) - b )*h^2;

   L1_u(count) = sum( abs( u(:) - b ) )*h^2;

   % computing maximum pointwise variation from previous u:
   
   variation_sup(count) = max( abs( u(:) - u_previous(:) ) );
   
   % computing L1-norm of solution variation from previous u:
   
   variation_L1(count) = sum( abs( u(:) - u_previous(:) ) )*h^2;
   
   % current solution values become previous solution values
   % for the next iteration cycle:
   
   u_previous = u;
   
   % printing out some statistics before going to next cycle:
   
   disp((t(count)-t0)/(tF-t0))
   disp([variation_L1(count), variation_sup(count), far_field_value(count)]) 
   
end
 
total_time_spent = sum( time_per_cycle )

figure(10)
plot(t,variation_L1,'-r')
xlabel('time')
ylabel('variation (L1 norm)')
title('L1 norm of last variation')

figure(20)
plot(t,variation_sup,'-r')
xlabel('time')
ylabel('variation (supnorm)')
title('supnorm of last variation')

figure(30)
plot(t,far_field_value,'-r')
xlabel('time')
ylabel('far field value')
title('far field value')

[X,Y] = meshgrid(x,y);

figure(50)
mesh(X,Y,u)
colormap(jet)
shading interp
xlabel('y')
ylabel('x')

%
% ----------------------------------------------------------------
%

disp(' ')
disp('*******************************************************************')
disp(' ')
disp('time per iteration cycle:')
disp(time_per_cycle')
disp(' ')
disp('*******************************************************************')
disp(' ')
disp('SUCCESSFUL EXECUTION')
disp(' ')
disp('*******************************************************************')

