function [ freq,Tr,Wr,lonData,latData,timeData ] = WeatherMapping( T,W,wave )
% %Author: Mozhou Gao
% %University of Calgary 
% %MGIS student in Department of Geography
% %Input: The wind, Temperatur, and precipitation threshold
% %Output: The 3D data box 


%% Data Loading
%  make sure data and the script are in the same folder
myFolder = 'C:\Users\mozhou\Desktop\1124Update';

if ~isdir(myFolder)% check the existence of folder 
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
%Locate all the .nc file in directory
TempPattern = fullfile(myFolder, '15YearWT.nc'); 
TempFiles = dir(TempPattern);
n = length(TempFiles);
PrepPattern = fullfile(myFolder, '15Yearprecip.nc');
PrepFiles = dir(PrepPattern);

TempFileName = TempFiles(1).name;
PrepFileName = PrepFiles(1).name;
latData = ncread(TempFileName,'latitude');
lonData = ncread(TempFileName,'longitude');
timeData = ncread(TempFileName,'time');
%% set the dimension of Data and preprocessing data 
la = length(latData);
lo = length(lonData);
m  = length(timeData);
Tr = zeros(m,1);
Wr = zeros(m,1);
Temp = ncread(TempFileName,'t2m');
Temp = Temp -273.15;%convert all the temperature to celcius degree  
uwnd = ncread(TempFileName, 'u10');
vwnd = ncread(TempFileName, 'v10');
Speed = sqrt(uwnd.^2 + vwnd.^2); %calculatet the net wind speed
Prep = ncread (PrepFileName,'tp');
%% Check the threshold of three parameters
boolT = zeros(lo,la,m);
boolW = zeros(lo,la,m);
boolP = zeros(lo,la,m);
for k = 1:1:m
    T1 = T-wave;
    T2 = T+wave;
    rt = (T2-T1).*rand(1,1)+T1;
    W1 = W-wave;
    W2 = W+wave;
    rw = (W2-W1).*rand(1,1)+W1;
    for i = 1:1:lo
        for j = 1:1:la
            if Speed(i,j,k)<=rw
                boolW(i,j,k) = 1;
            end
            if Temp(i,j,k)>=rt
                boolT(i,j,k) = 1;
            end
            if Prep(i,j,k) <=0
                boolP (i,j,k) =1;
            end
        end
    end
    Tr(k) = rt;
    Wr(k) = rw;
end
newbool = boolW + boolT + boolP;
C = (newbool == 3);
  
%% Output
freq =C;

end

