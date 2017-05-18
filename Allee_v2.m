%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Gillespie Allee v2   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Walter Ullon 6/17/16
% This version generates a predetermined number of time series given by the
% index of the 1st 'for' loop. Each time series is saved to a CSV file.
% Whenever the population reaches extinction, the prefix "EX" is added to
% the file name in order to differentiate it from those in which no
% extinction takes place.
% Code was adapted from Allee.M. All unnecessary pieces of code 
% were deleted to streamline the simulations.

for i=1:50
SEED = i;
stream0=RandStream('mt19937ar','Seed', i);
RandStream.setGlobalStream(stream0)

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
        
        break
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Write to csv file        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     M = [t_new; X_new]';                          % matrix of time vs. pop (for CSV)
%     if X <= 0
% %         fileName1 = ['Ex_Allee' num2str(SEED) '-N120.csv'];
%         fileName = ['Allee' num2str(SEED) '-N120.csv'];
%         csvwrite(fileName, M);                        % write to CSV
% %         csvwrite(fileName1, M);                        % write to CSV
%     else
%         fileName = ['Allee' num2str(SEED) '-N120.csv'];
%         csvwrite(fileName, M);                        % write to CSV
%     end
%   
    
    %disp(SEED)
    disp(t_extinct)
  
    clear all
end

   