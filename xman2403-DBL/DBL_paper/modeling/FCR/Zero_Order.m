
% Inputs are:   DO concentration (mg/L)
%               Vertical position (mm)
%               Time of each profile (Serial date)
%               Diffusion coefficient (m2/s)
% Set up variables
 load 'epI_sed_0_20210712.mat'; % 16/02/2021
 DO_sed(:,40) =0;
 load 'time_swi_epi_20210216_2pm.mat'; 
 load 'depth_lander_201030.mat';
kmin=100;
kstep=50;
kmax=800;
t_total=time_swi;
k_vals = kmin:kstep:kmax;       % Range of first-order rate constants
k_vals = k_vals/86400;% Convert from d^-1 to s^-1
DO = DO_sed';
DO_sed(DO_sed<=0.096) = 0;
BO = DO_sed';
% DO = DO_sed';
[row_total,column] = size(DO);        % Get number of rows and columns of DO
 Z=position_sed;
 %Ds=1.213e-09; % CCR
 Ds=1.236e-09;  % FCR epi
 
  profile_times = zeros(1,length(t_total)); %length(t)=174
  profile_times(1) = 0;
    for a = 2:length(t_total)     
         profile_times(a) = (t_total(a)-t_total(a-1));  % times in seconds
    end
         k_profile = zeros (1,length(profile_times));
         rls_best_profile = zeros (1,length(profile_times));
         
   %read row of data for each time profile
 row = zeros (1,length(profile_times));
 for m = 1:column
 for n =1:row_total
     if DO(n,m)> 0
         row (1,m) = row (1,m)+1;
     end
 end
     row (1,m) = row (1,m)+1;
 end
 
 row_all= max(row(1,:));
 DO_calc_best = zeros (row_all,column);
 
  bot = zeros (1,length(profile_times));
 for m = 1:column
 for n =1:row_total
     if BO(n,m)> 0
         bot (1,m) = bot (1,m)+1;
     end
 end
     bot (1,m) = bot (1,m)+1;
 end
 
 
 for j = 1:length(profile_times)-1
  DO_calc_best_each = zeros(row(1,j),1);   % Initialize DO_calc_best variable (stores best DO profiles)
  rls_best = 1e30;                % Initialize R2 variable 
  dt = 1;                         % Time step (s) for caclulated profiles
  dz = 0.1/1000;  % dz = 0.1mm; converts to m
  t=time_swi(j:j+1);
  times = time_swi(j+1)-time_swi(j);
  steps = times/dt;               % Calculate number of time steps to simulate
  profile_dum = zeros (40,1);     % the modeled profile for RMSE calculation 
  DO_dum = zeros (40,1);

%Let row be the bottom boundary of DO_sed (Z)
   DO_calc_1 = zeros(row(1,j),steps); 
   DO_calc_1(:,1) = DO(1:row(1,j),j);
    %Transverse the bottom and upper boundary
   
% Add upper boundary condition from field data and interpolate to fill in
% the boundary conditions on the additional simulated profiles
      if j == 1
        i1 = time_swi(j)/dt+1;
      else
        i1 = time_swi(j)/dt;    % Step # of i-1 profile
      end
    i2 = time_swi(j+1)/dt;      % Step # of i profile
    c1 = DO(1,j);                 % DO conc of i-1 profile
    c2 = DO(1,j+1);                   % DO conc of i profile
    dummy = nan(1,(i2-i1-1));       
      for l = 1:length(dummy)       % Interpolate Cswi between i-1 and i profiles, store in dummy
        dummy(l) = c1+(c2-c1)/(i2-i1)*l;
      end
    DO_calc_1(1,2:(i2-i1)) = dummy;   % Fill in DO_calc_1 with interpolated Cswi values
    DO_calc_1(1,i2-i1+1) = DO(1,j+1);      % Add next measured Cswi to DO_calc_1 variable
 
 %Set bottom boundary condition
       i = 1;
   for i = 1:length(DO_calc_1)
    if row(1,j)>= 2
      DO_calc_1(row(1,j),i) = DO_calc_1((row(1,j)-1),i) - 0.4 ;
    end
   
    if bot(1,j+1) >= row(1,j)
       if DO_calc_1(row(1,j),i)< 0    % whether to reset to zero
          DO_calc_1(row(1,j),i) = 0;
       end  
    end 
   end

% The Meat & Potatoes of the code!
   for a = 1:length(k_vals) % Loop over each value of k
    
    k = k_vals(a); % Set First-order rate constant
    C = DO_calc_1; % Start with fresh DO matrix
    rls_profile = zeros(1,1);
       for i = 2:steps          % Loop over every dt step
            penetration = row (1,j)-1;
           for n = 2:(row(1,j)-1)  % Loop over each depth position
         %for n = 2:18;           % Loop over each depth position
         C(n,i) = (Ds*dt/(dz^2))*(C(n+1,i-1)-2*C(n,i-1)+C(n-1,i-1))-k*dt+C(n,i-1);
                 if (C(n-1,i)>=0 && C(n,i)<=0)
                      penetration = n-1;
                 end
               
                 if n<= bot(1,j+1)
                    if C(n,i) < 0 % Force to go to zero if results in negative concentration
                       C(n,i) = 0;
                    end
                 end
            end     % end depth loop (n)
             C(penetration+1,i) = C(penetration,i) - 0.4;
             if n <= bot(1,j+1)
                if C(n,i) < 0 % Force to go to zero if results in negative concentration
                      C(n,i) = 0;
                end
             end
        end         % end profile loop (i)
    
     % calculate RMSE for the whole profile % 24/02/2021
            pos = profile_times(j+1);  % /2 in old
            
            if bot(j)>=bot(j+1)
                profile_dum (1:bot(j+1),1) = C(1:bot(j+1),pos); % should be the last column of C 
            else
                profile_dum (1:bot(j),1) = C(1:bot(j),pos); 
            end
             DO_dum (:,1) = DO(:,j+1);  % field DO profile 24/02/2021
             rls=sqrt((1/bot(j))*sum((DO_dum(:)-profile_dum(:,1)).^2));
             rls(isnan(rls)) = 0;  % Remove NaNs that arise if dividing by zero
             rls(isinf(rls)) = 0;  % Remove Infs that arise if dividing by zero
    
    % Compare R2 to k_best
    if rls < rls_best
        rls_best_profile(1,j) = rls;
        rls_best = rls;
        k_profile(1,j) = k;
        % Save C matrix if best DO concentrations
        DO_calc_best_each(:,1) = C(:,(profile_times(j+1)));
        DO_calc_best(1:row(1,j)-1,j+1) = DO_calc_best_each (1:row(1,j)-1,1);
    end    

   end % end k loop (a)
  k_best = sum(k_profile(1:length(profile_times)))/length(k_profile);
  rls_best = sum(rls_best_profile(1:length(profile_times)))/length(rls_best_profile);
 end

DO_calc_best (1:row(1,1),1) = DO(1:row(1,1),1);
k_best=k_best*86400;    % Convert back to d^-1;
