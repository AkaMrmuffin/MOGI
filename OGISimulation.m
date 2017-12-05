function [Tr,Wr,daily_reparied_wells,wells_missed_ripeness,wells_missed_weather...
    ,prop_cells_unavail,daily_missed_visits,daily_emissions,daily_leaks,...
    Total_emissions,daily_repaired_leaks,outstanding_wells,freq,dailycell_Emissions,daily_new_leaks,cnt] = OGISimulation(T,W,Wells_Day,Min_Int,sampnum,wave)
%Weather Mapping
[ freq,Tr,Wr,lonData,latData,timeData ] = WeatherMapping( T,W,wave );
t = length(timeData);
lat = latData;
long = lonData;
m = length(lat);
n = length(long);
%% load the test well data
num = xlsread('Wells.xlsx'); 
rnum = datasample(num,sampnum); 
wlat = rnum(:,1);
wlon = rnum(:,2);   
wlon360 = mod (wlon, 360.0);        
%% count number of wells in grid
cnt = hist3([wlon360,wlat], {240:0.125:250  49:0.125:60});

%% Leak Initiation
Init_leaks  = Initializ( cnt,n,m );
%% FWAQS
MAmount = xlsread('FWAQS_distribution.xls'); %read the Methane Emission Size data 
Msize = MAmount(:,4); 
SimMsize = zeros(n,m); 
%random sample the number of Emission sizes from FWAQS data for each generated leaks
for i =1:1:n
    for j = 1:1:m
        if Init_leaks(i,j) ~=0
           SimMsize(i,j) = sum(datasample(Msize,Init_leaks(i,j)).*86400);
        end 
    end 
end 

%% Define Some Important variables  
t_last_LDAR = zeros(n, m);% record the time when there was last LDAR of each cell 
outstanding_wells = zeros(n, m); %record number of wells we couldn't complete in each grid cell
dailycell_Emissions = zeros(n,m,t);
Total_emissions = zeros(t, 1); %The total emssion size for alberta for every day
daily_emissions = zeros(t, 1);
daily_leaks = zeros(t,1);
daily_new_leaks = zeros(t,1);
daily_repaired_leaks = zeros(t,1);
prop_cells_unavail = zeros(t,1);% proportion of province where weather is not ok
%vector record The wells can't fix due to specific reasons
wells_missed_ripeness = zeros (t, 1); 
wells_missed_weather = zeros (t, 1);
% evaluate regulatory compliance by tracking when things aren't getting done to spec
noncompliance = zeros (n, m, t);     % full 3d array
daily_missed_visits = zeros (t, 1);     % just the timeseries
%% new parameters that should be defined outside of this function and input as args
num_wells_per_day = Wells_Day;  % this is the number of wells that can be surveyed per day
minimum_interval = Min_Int;% this is the minimum allowable survey interval
max_LDAR_gap = 120;   
max_LDAR_gap_matrix = max_LDAR_gap.*ones(n,m);
Pnew = 0.00133; % probability of having new leak for each well
% num_wells_per_day = number_LDAR_crews * (working_hours / time_to_survey_one_well)
% cost = number_LDAR_crews * cost_per_day_per_LDAR_crew - when everyone is working flat out
startsim_time = zeros(t,1);
daily_reparied_wells = zeros(t,1);
for k = 1:1:t
    Nleaksbool = ( cnt > 0);
    Pnewleaks = Nleaksbool.*Pnew;% create a probability matrix 
    Newleaks = binornd(cnt,Pnewleaks); %generate a new leak for each day of the year
    NMsizeold = zeros(n,m);
    NMsizenew = zeros(n,m);
    %Emssion SIze Calculateion.
    for i = 1:1:n
        for j =1:1:m
               NMsizeold(i,j) = sum(datasample(Msize,Init_leaks(i,j)).*86400);
               NMsizenew(i,j) = sum(datasample(Msize,Newleaks(i,j)).*86400);
         end
    end
    SimMsize = SimMsize+NMsizeold +NMsizenew;
    daily_emissions(k) = sum(sum(NMsizenew + NMsizeold));
    dailycell_Emissions(:,:,k) = NMsizenew + NMsizeold;
    daily_leaks(k) = sum(sum(Init_leaks));
    daily_new_leaks(k) = sum(sum(Pnewleaks));
    Total_emissions (k) = sum (sum(SimMsize));
    Init_leaks = Init_leaks + Newleaks; %add the newleaks to the number of leaks
    prop_cells_unavail (k) = 1 - (sum (sum (freq(:,:,k))) / (n * m));
    
    Morewells = true;
    num_wells_remaining_today = num_wells_per_day;
    while Morewells
        t_since_last_LDAR = k - t_last_LDAR;
        % reshape to linear array and sort by time to get an array of all the most overdue cells, sorted descending
        t_since_last_LDAR_linear = t_since_last_LDAR (:);
        % Shuffling
        [t_since_last_LDAR_sorted, indices] = sort (t_since_last_LDAR_linear, 'descend');
        [ t_since_last_LDAR_sorted, indices] = ShuffleIndex( t_since_last_LDAR_sorted, indices );
       
        %%%%spin%%%%
        for i = 1:length(t_since_last_LDAR_sorted)
            % first check to see if this cell has been long enough since last visit
             if t_since_last_LDAR_sorted (i) < minimum_interval
                 Morewells = false;
                 wells_missed_ripeness (k) = num_wells_remaining_today;
                 break;             % break out of this loop, we don't need to check the rest as
                                    % they will all have less time
             else
                [n_target, m_target] = ind2sub (size (freq(:,:,k)), indices (i));
                if freq (n_target, m_target, k) == 1 && cnt (n_target, m_target) > 0
                   if outstanding_wells (n_target, m_target) > 0
                       num_wells = outstanding_wells (n_target, m_target);%we have some oustanding wells  
                       daily_reparied_wells(k) =  num_wells;
                   else
                       num_wells = cnt(n_target,m_target); %all wells
                       if num_wells > num_wells_per_day
                            daily_reparied_wells(k) =  num_wells_per_day;
                       else
                            daily_reparied_wells(k) = num_wells;
                       end
                   end
                   if num_wells_remaining_today > num_wells
                       num_wells_remaining_today = num_wells_remaining_today - num_wells;
                       t_last_LDAR (n_target, m_target) = k;       % record to last LDAR date
                       outstanding_wells (n_target, m_target) = 0; % we have completed all wells here
                       
                   else
                      
                       outstanding_wells (n_target, m_target) = num_wells - num_wells_remaining_today;
                       num_wells_remaining_today = 0;
                       Morewells = false;
                   end
                   % mitigate some emissions, repair those leaks. First calculate the fraction we completed
                   
                   Frac_of_wells_LDar = 1 - (outstanding_wells(n_target,m_target)/num_wells);                 
                   Init_leaks(n_target,m_target) = Init_leaks(n_target,m_target) - round((Init_leaks(n_target,m_target)*Frac_of_wells_LDar));
                   SimMsize(n_target, m_target) = SimMsize(n_target, m_target) - (SimMsize(n_target, m_target)*Frac_of_wells_LDar);
                   daily_repaired_leaks(k) = round((Init_leaks(n_target,m_target)*Frac_of_wells_LDar));
                   
                elseif freq (n_target, m_target, k) == 0 && cnt (n_target, m_target) > 0
                   well_ava = num_wells_per_day;
                   if outstanding_wells (n_target, m_target) == 0 ||  outstanding_wells (n_target, m_target) > well_ava
                       wells_missed_weather(k) = well_ava;
                   else    
                       wells_missed_weather(k) = outstanding_wells (n_target, m_target);
                   end
                else 
                  continue
                end
             end
        end
        
        Morewells = false; %flag, delete it after all  
    end
    
    %figure out how compliant you are
    % check to assign the starting date of the simulation (the simulation cannot be used without spinup, you must LDAR every cell to reset the clock per se)
    % and allow the t_last_LDAR to be a real number (not 0). In reality the t_last_LDAR will be in the year before the start of the simulation, but we don't
    % have that data, so you can't use these data.
    time_difference_cell = k - t_last_LDAR;
    noncompliance(:,:,k) = time_difference_cell > max_LDAR_gap_matrix;
    daily_missed_visits (k) = sum(sum(noncompliance(:,:,k)));
    zeroC = find (cnt > 0 & t_last_LDAR ~= 0);
    startsim_time(k) = length(zeroC);
    %sunset
    fprintf (strcat ('completed day: ', num2str(k), '\n'))
end

end

