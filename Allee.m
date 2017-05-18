clear
format long

stream0=RandStream('mt19937ar','Seed', 23);
RandStream.setGlobalStream(stream0)

fileName = 'Allee15_120.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Gillespie Allee    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% * This version captures the extinction time, plus it counts the number of
%   times the infected population has gone extinct.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                             
%                                                                         %
%                 X_dot = x*(x-x1)*(x-x2)                                 %
%                 X_dot = (-sigma/6)*x^3 + (lambda/2)*x^2 -mu*x           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%
n = 10^7;                           % number of iterations

N = 120;                            % Total Population (dN/dt = 0)
mu = 0.2;                           % Low density death rate.
siggma = 3.0;                       % Overcrowded negative growth rate (-).
lammda = 1.425;                     % Optimal density growth rate.
K = 100;                            % Carrying capacity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Steady States      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_1 = (3*lammda + sqrt(9*lammda^2 - 24*siggma*mu))/(2*siggma)*K;  % stable eq.
X_2 = (3*lammda - sqrt(9*lammda^2 - 24*siggma*mu))/(2*siggma)*K;  % unstable eq.

% X_dot1 = (-siggma/6)*(X_1)^3 + (lammda/2)*(X_1)^2 -mu*(X_1);
% X_dot2 = (-siggma/6)*(X_2)^3 + (lammda/2)*(X_2)^2 -mu*(X_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Initial conditions  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;                                           %%%  Time  %%%

X = N;

t_new(1) = 0;                                    %% Time array (even steps).                        
X_new(1) = X;                                    %% Population array (time course) %%                           

events = zeros(3,1);                                  %% Events array.

jump = 1;                           %% Designated time step.
j = 1;                              %% Time step index. 

t_extinct = 0;                      %% Captures the extinction time.
e_counter = 0;                      %% Counts the number of times gone extinct.

for i=1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw random numbers %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    r1=rand(1);
    r2=rand(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Events        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    events(1) = mu*X;                               % Low Density death      -------> 
    events(2) = lammda*X*(X-1)/(2*K);               % Optimal density growth    ---->  
    events(3) = (siggma*X*(X-1)*(X-2))/(6*K^2);     % Death due to Overcrowding ---->  
    
    a0 = sum(events);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Time step       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau = (1/a0)*log(1/r1);
    t = t + tau;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Actuating Gaussian event      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gauss = 0;
    atemp = 0;
    for k = 1:3
        atemp = atemp+events(k);
        if r2*a0 < atemp,
            gauss = k;
            break;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Update population      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if gauss == 1
        X = X - 1;
    elseif gauss == 2
        X = X + 1;
    else
        X = X - 1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% * Saves time values in equal-sized steps (jump).   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t > jump
            t_new(j+1) = floor(t);
            X_new(j+1) = X;
            
%%%%%%%%%%%%%%%%%%%%
%%%%%% * PULSE.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                        
%          if t_new(j+1) == 9901
%              p = 55;                           % Perturbation size
%              X = X - p;
%              
%              X_new(j+1) = X;
%          end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            jump = jump + 1;
            j = j + 1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Capture extinction     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if X<=0
        X_new(j+1) = X;
        t_new(j+1) = floor(t);
        
        e_counter = e_counter + 1;
        t_extinct = t_extinct + floor(t);
        disp(t_extinct)
        break
    end
    
end


M = [t_new; X_new]';                          % matrix of time vs. pop (for CSV)
% csvwrite(fileName, M);                        % write to CSV


X_new = X_new';


    hold on
    plot(t_new, X_new);
    plot([t_new(1), t_new(end)],[K, K], 'g', 'Linewidth', 1);
    plot([t_new(1), t_new(end)],[X_1, X_1], 'm', 'Linewidth', 1);
    plot([t_new(1), t_new(end)],[X_2, X_2], 'm--', 'Linewidth', 1);
%     plot(9901,104,'ro','Markersize', 12, 'Linewidth',1.5);
    axis([0,3.5*10^4,0,180])
%     title('Gillespie with Allee effect');
    xlabel('Time')                              % t-axis label
    ylabel('Population')                        % N-axis label
%     legend('population', 'carrying capacity', 'stable equilibrium', 'unstable equilibrium','high resilience')
%     legend('population', 'carrying capacity')
