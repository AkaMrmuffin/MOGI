Trange = -3;
Wrange = 10;
Lt = length(Trange);
Lw = length(Wrange);
Time = ncread('15YearWT.nc','time');
lat = ncread('15YearWT.nc','latitude');
long = ncread('15YearWT.nc','longitude');
lo = length(long);
la = length(lat);
day = [ 365 365 365 366 365 365 365 366 365 365 365 366 365 365 365 366];
lengthTime = length(Time);
T_Year = [];
T_T = [];
T_W = [];
T_lat = [];
T_long = [];
T_OGIdays = [];
T_MWGs =[];
T_AnnualEmission = [];
TimeID = [];
TT = [];
TW = [];
T_daily_reparied_wells = [];
T_wells_missed_ripeness = [];
T_wells_missed_weather = [];
T_prop_cells_unavail = [];
T_daily_missed_visits = [];
T_daily_emissions = [];
T_daily_leaks =[];
T_Total_emissions = [];
T_daily_repaired_leaks = [];
for i = 1:1:Lt
    for j = 1:1:Lw
        T = Trange(i);
        W = Wrange(j);
        index1 = 1;
        index2 = 365;
        %%  Main Model
       [Tr,Wr,daily_reparied_wells,wells_missed_ripeness,wells_missed_weather...
    ,prop_cells_unavail,daily_missed_visits,daily_emissions,daily_leaks,...
    Total_emissions,daily_repaired_leaks,outstanding_wells,freq,dailycell_Emissions,daily_new_leaks] = OGISimulation(T,W,10,60,5000,0.5);
      %Spatial Based
       OGIday = zeros(lo,la,length(day));
       MWGsday = zeros(lo,la,length(day));
       Annual_Emission = zeros(lo,la,length(day));
       temp_Tr = ones(lo,la,length(day));
       temp_Wr = ones(lo,la,length(day));
       temp_year = ones(lo,la,length(day));
       temp_lat = zeros(lo,la,length(day));
       temp_long = zeros(lo,la,length(day));
       for k = 1:1:length(day)
                YearOGI = freq(:,:,index1:index2);
                Annual_Emission(:,:,k) = sum(dailycell_Emissions(:,:,index1:index2),3);
                OGIday(:,:,k) = sum(YearOGI,3);
                temp_Tr(:,:,k) = mean(Tr(index1:index2));
                temp_Wr(:,:,k) = mean(Wr(index1:index2));
                [X,Y] = meshgrid(long,lat);
                temp_lat(:,:,k) = Y';
                temp_long(:,:,k) = X';
                temp_year(:,:,k) =temp_year(:,:,k).* k;
                for l = 1:1:lo
                    for m = 1:1:la
                        %find the index of evary OGI working day 
                         OGIdayIndex = find(YearOGI(l,m,:));
                         L =length(OGIdayIndex);%check the number of OGI days for each cell
                         if L == 0 %if OGI days = 0 
                            MWGsday(l,m,k)=365 ;
                         elseif L ==1 %if OGI days = 1
                             tempM = YearOGI(l,m,:); 
                            MWGsday(l,m,k)=365 - OGIdayIndex; %Max Winter gap = 366
                         else %otherwise
                             Diff = diff(OGIdayIndex);% calculate the differences between every OGI days index 
                             FD = OGIdayIndex(1);
                             LD = 365 - OGIdayIndex(L);
                             MaxGap = max(Diff);%Find the maximum difference
                             TempList = [FD LD MaxGap];
                             MWGsday(l,m,k) = max(TempList);
                         end 
                    end
                end
                index1 = index1 + day(k);
                index2 = index2 + day(k);                  
       end
       T_Year = vertcat(T_Year,temp_year(:));
       T_T = vertcat(T_T,temp_Tr(:));
       T_W = vertcat(T_W,temp_Wr(:));
       T_lat = vertcat(T_lat,temp_lat(:));
       T_long = vertcat(T_long,temp_long(:));
       T_OGIdays = vertcat(T_OGIdays,OGIday(:));
       T_MWGs =vertcat(T_MWGs,MWGsday(:));
       T_AnnualEmission = vertcat(T_AnnualEmission,Annual_Emission(:));
       %Time Series Based 
        TT= vertcat(TT, Tr);
        TW= vertcat(TW, Wr);
        T_daily_reparied_wells =vertcat(T_daily_reparied_wells, daily_reparied_wells);
        T_wells_missed_ripeness = vertcat(T_wells_missed_ripeness,wells_missed_ripeness);
        T_wells_missed_weather = vertcat(T_wells_missed_weather,wells_missed_weather);
        T_prop_cells_unavail = vertcat(T_prop_cells_unavail,prop_cells_unavail);
        T_daily_missed_visits = vertcat(T_daily_missed_visits,daily_missed_visits);
        T_daily_emissions = vertcat(T_daily_emissions,daily_emissions);
        T_daily_leaks = vertcat(T_daily_leaks,daily_leaks);
        T_Total_emissions = vertcat(T_Total_emissions,Total_emissions);
        T_daily_repaired_leaks = vertcat(T_daily_repaired_leaks,daily_repaired_leaks);
        time = 1:lengthTime;
        TimeID = horzcat(TimeID, time);

       fprintf (strcat ('completed Simulation: ', num2str(i+j), '\n'))
    end
end
%Output The Data
Output = [T_Year, T_T ,T_W ,T_lat,T_long,T_OGIdays,T_MWGs,T_AnnualEmission];
Table = array2table(Output,'VariableNames',{'Year','T_thresh_day','W_thresh_day','latitude'...,
                                            'longitude','OGIdays','MaxWinterGap','AnnualEmission'});
                                        
 writetable(Table,'Output.csv');
 
Output2 = [TimeID', TT, TW ,T_daily_reparied_wells,T_wells_missed_ripeness,...
    T_wells_missed_weather,T_prop_cells_unavail,T_daily_missed_visits,...
    T_daily_emissions,T_daily_leaks,T_Total_emissions,T_daily_repaired_leaks];
Table2 = array2table(Output2,'VariableNames',{'TimeID','T_thresh_day','W_thresh_day','daily_reparied_wells',...
                                            'wells_missed_ripeness','wells_missed_weather','prop_cells_unavail'...
                                            'daily_missed_visits','daily_emissions','daily_leaks','Total_emissions','daily_repaired_leaks'});
writetable(Table2,'Output2.csv');          



