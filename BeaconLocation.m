function [solutions, recalculate] = BeaconLocation(freq)


clc; clear all; close all;


sensorArray = audioDeviceReader("Device","IN 1-4 (2- BEHRINGER UMC 404HD 192k)", ...
    "Driver","DirectSound","NumChannels",4,"SampleRate",192000,"BitDepth",'16-bit integer')

% audioFromDevice=sensorArray() % gets one frame of data from device.
                                % samples per frame can be specified, is
                                % the "buffer" which can be defined and is
                                % the device record start time latency.
                                
setup(sensorArray)

fileWriter = dsp.AudioFileWriter('hydrophones.wav','FileFormat','WAV',"SampleRate",192000);

disp("Record start")

tic

while toc<2.5
    
    acquiredAudio=sensorArray();
    
    fileWriter(acquiredAudio);
    
end

disp("Record end")

release(sensorArray)

release(fileWriter)

[data, Fs]=audioread('hydrophones.wav');

t = 0:1/Fs:(length(data)-1)/Fs;

sensor_A = data(:,1);
sensor_B = data(:,2);
sensor_C = data(:,3);
sensor_D = data(:,4);

fs = 192000 % Set sampling rate to variable "Fs"

tA = linspace(0,length(sensor_A)/fs,length(sensor_A));
tB = linspace(0,length(sensor_B)/fs,length(sensor_B));
tC = linspace(0,length(sensor_C)/fs,length(sensor_C));
tD = linspace(0,length(sensor_D)/fs,length(sensor_D));

F_sensor_A = fft(sensor_A); % FFT on "sensor_A"
F_sensor_B = fft(sensor_B);
F_sensor_C = fft(sensor_C);
F_sensor_D = fft(sensor_D);
    
fvec_sensor_A = linspace(-fs/2,fs/2,length(sensor_A)); % frequency vector
fvec_sensor_B = linspace(-fs/2,fs/2,length(sensor_B));
fvec_sensor_C = linspace(-fs/2,fs/2,length(sensor_C));
fvec_sensor_D = linspace(-fs/2,fs/2,length(sensor_D));

% Filter audio, Remove errant peaks, create plots. 
% Filtered & cleaned audio from each sensor A is: 
% "clean_filtered_A"
% "clean_filtered_B"
% "clean_filtered_C"
% "clean_filtered_D"

% Filter for sensor A

[b,a] = butter(5,[freq*1000-750, freq*1000+750]/(Fs/2),'bandpass');
    filtered_A = filter(b,a,sensor_A);
    FFT_filtered_A = fftshift(abs(fft(filtered_A)));
    clean_filtered_A = filtered_A;
    
    thresh = .005;
    count = 0;
    
    W = 3000;                               % Window width
    for i = W+1:length(clean_filtered_A)-W
%     for i = 175700:175740                 Manually specified test window
        if abs(clean_filtered_A(i))>thresh
            count=sum(abs(clean_filtered_A(i-W:i+W))>thresh);
            if count<0.2*W
                out2(i)=0;
            end
        end
    end
     
% Filter for sensor B

[b,a] = butter(5,[freq*1000-750, freq*1000+750]/(Fs/2),'bandpass');
    filtered_B = filter(b,a,sensor_B);
    FFT_filtered_B = fftshift(abs(fft(filtered_B)));
    clean_filtered_B = filtered_B;
    
    thresh = .005;
    count = 0;
    
    
    for i = W+1:length(clean_filtered_B)-W
%     for i = 175700:175740
        if abs(clean_filtered_B(i))>thresh
            count=sum(abs(clean_filtered_B(i-W:i+W))>thresh);
            if count<0.2*W
                out2(i)=0;
            end
        end
    end
  
% Filter & Plots for sensor C

[b,a] = butter(5,[freq*1000-750, freq*1000+750]/(Fs/2),'bandpass');
    filtered_C = filter(b,a,sensor_C);
    FFT_filtered_C = fftshift(abs(fft(filtered_C)));
    clean_filtered_C = filtered_C;
    
    thresh = .005;
    count = 0;
    
    
    for i = W+1:length(clean_filtered_C)-W
%     for i = 175700:175740
        if abs(clean_filtered_C(i))>thresh
            count=sum(abs(clean_filtered_C(i-W:i+W))>thresh);
            if count<0.2*W
                out2(i)=0;
            end
        end
    end

% Filter & Plots for sensor D

[b,a] = butter(5,[freq*1000-750, freq*1000+750]/(Fs/2),'bandpass');
    filtered_D = filter(b,a,sensor_D);
    FFT_filtered_D = fftshift(abs(fft(filtered_D)));
    clean_filtered_D = filtered_D;
    
    thresh = .005;
    count = 0;
    
    
    for i = W+1:length(clean_filtered_B)-W
%     for i = 175700:175740
        if abs(clean_filtered_D(i))>thresh
            count=sum(abs(clean_filtered_B(i-W:i+W))>thresh);
            if count<0.2*W
                out2(i)=0;
            end
        end
    end
    


% TDOA algorithm that uses filtered audio to determine timestamp of desired signal arrival,
% then trims unfiltered audio array around that timestamp. The arrival of
% the unfiltered audio to each hydrophone is used to determine TDOA.

%   Sensor Pair A & B:

tempA = find(clean_filtered_A >.001, 1); % Find the array index of clean_filtered_A
                                         % where amplitude first begins to climb.

TDOA_leftedge  = 210;
TDOA_rightedge = 521;

TDOAsensorA = abs(sensor_A(tempA-TDOA_leftedge:tempA+TDOA_rightedge)); % New trimmed unfiltered audio array

% Repeat for sensor_B

TDOAsensorB = abs(sensor_B(tempA-TDOA_leftedge:tempA+TDOA_rightedge));

DiffThreshA = .15*(max(TDOAsensorA));
DiffThreshB = .15*(max(TDOAsensorB));

TDOAindex_A = find(TDOAsensorA>DiffThreshA,1);
TDOAindex_B = find(TDOAsensorB>DiffThreshB,1);

TDOA_time_A = t_TDOApair(TDOAindex_A);
TDOA_time_B = t_TDOApair(TDOAindex_B);

TDOA_time_AB = TDOA_time_A-TDOA_time_B;

tau_AB = round((TDOA_time_A-TDOA_time_B)*fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sensor pair A & C:

% Repeat for sensor_C

TDOAsensorC = abs(sensor_C(tempA-TDOA_leftedge:tempA+TDOA_rightedge));

DiffThreshA = .15*(max(TDOAsensorA));
DiffThreshC = .15*(max(TDOAsensorC));

TDOAindex_A = find(TDOAsensorA>DiffThreshA,1);
TDOAindex_C = find(TDOAsensorC>DiffThreshC,1);

TDOA_time_A = t_TDOApair(TDOAindex_A);
TDOA_time_C = t_TDOApair(TDOAindex_C);

TDOA_time_AC = TDOA_time_A-TDOA_time_C;

tau_AC = round((TDOA_time_A-TDOA_time_C)*fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sensor pair A & D

% Repeat for sensor_D

TDOAsensorD = abs(sensor_D(tempA-TDOA_leftedge:tempA+TDOA_rightedge));

DiffThreshA = .15*(max(TDOAsensorA));
DiffThreshD = .15*(max(TDOAsensorD));

TDOAindex_A = find(TDOAsensorA>DiffThreshA,1);
TDOAindex_D = find(TDOAsensorD>DiffThreshD,1);

TDOA_time_A = t_TDOApair(TDOAindex_A);
TDOA_time_D = t_TDOApair(TDOAindex_D);

TDOA_time_AD = TDOA_time_A-TDOA_time_D;

tau_AD = round((TDOA_time_A-TDOA_time_D)*fs);

% End unfiltered audio time-lag component

% Begin fine tuning TDOA values with Cross Correlation

sensor_A = data(:,1);
sensor_B = data(:,2);
sensor_C = data(:,3);
sensor_D = data(:,4);

w=clean_filtered_A; 
x=clean_filtered_B;
y=clean_filtered_C;
z=clean_filtered_D;

lA = length(w);
lB = length(x);
lC = length(y);
lD = length(z);

samples1 = 1:min(lA,lB);
samples2 = 1:min(lA,lC);
samples3 = 1:min(lA,lD);


[C1, lag1] = xcorr(w(samples1), x(samples1), 155);
[C2, lag2] = xcorr(w(samples2), y(samples2), 155);
[C3, lag3] = xcorr(w(samples3), z(samples3), 135);

%normalizing XCorr
C1 = C1/max(C1);
C2 = C2/max(C2);
C3 = C3/max(C3);

WINDOW = ceil(fs/(2*1000*freq)*1.5);

n1 = ceil(length(C1)/2);
n2 = ceil(length(C2)/2);
n3 = ceil(length(C3)/2);

    
C1sub2 = C1((tau_AB + n1 - WINDOW):(tau_AB + n1 + WINDOW));
idx_AB2 = find(C1sub2==max(C1sub2))+tau_AB+n1-WINDOW;

C2sub2 = C2((tau_AC + n2 - WINDOW):(tau_AC + n2 + WINDOW));
idx_AC2 = find(C2sub2==max(C2sub2))+tau_AC+n2-WINDOW;

C3sub2 = C3((tau_AD + n3 - WINDOW):(tau_AD + n3 + WINDOW));
idx_AD2 = find(C3sub2==max(C3sub2))+tau_AD+n3-WINDOW;


tau_AB2 = lag1(idx_AB2);
tau_AC2 = lag1(idx_AC2);
tau_AD2 = lag1(idx_AD2);

C1sub = C1((tau_AB2 + n1 - WINDOW):(tau_AB2 + n1 + WINDOW));
idx_AB = find(C1sub==max(C1sub))+tau_AB2+n1-WINDOW;

C2sub = C2((tau_AC2 + n2 - WINDOW):(tau_AC2 + n2 + WINDOW));
idx_AC = find(C2sub==max(C2sub))+tau_AC2+n2-WINDOW;

C3sub = C3((tau_AD2 + n3 - WINDOW):(tau_AD2 + n3 + WINDOW));
idx_AD = find(C3sub==max(C3sub))+tau_AD2+n3-WINDOW;


adjusted_tdoa_seconds_A_B = lag1(idx_AB)/fs
adjusted_tdoa_seconds_A_C = lag1(idx_AC)/fs
adjusted_tdoa_seconds_A_D = lag1(idx_AD)/fs


% Multilateration now needs to take the above three
% adjusted_tdoa_seconds_A_X" variables as inputs





%% Speed of sound in water
c = 1500;        % m/s
%% Location of the hydrophones
% Sensor A is at the origin. Do not change this.
xA = 0;
yA = 0;
zA = 0;

xB = -.838;
yB = -.813;
zB = .254;


xC = -.838;
yC = -.406;
zC = .813;

xD = -.838;
yD = -.406;
zD = .419;

syms x y z

eq1 = -round(TDOA_time_AB, 5) == (1/c) * ( sqrt( (x-xB)^2 + (y-yB)^2 + (z-zB)^2 ) - sqrt(x^2 + y^2 + z^2) );
eq2 = -round(TDOA_time_AC, 5) == (1/c) * ( sqrt( (x-xC)^2 + (y-yC)^2 + (z-zC)^2 ) - sqrt(x^2 + y^2 + z^2) );
eq3 = -round(TDOA_time_AD, 5) == (1/c) * ( sqrt( (x-xD)^2 + (y-yD)^2 + (z-zD)^2 ) - sqrt(x^2 + y^2 + z^2) );

sol = solve([eq1, eq2, eq3], [x, y, z]);

xSol = double(sol.x);
ySol = double(sol.y);
zSol = double(sol.z);

if length(xSol) == 1
    solutions = [xSol,ySol,zSol]
    recalculate = 0; 
elseif isempty(xSol) 
    solutions = [0,0,0,0,0,0]; %% Solution using orignal TDOA, no XCorr
    recalculate = 1; 
else    
    solutions = [xSol(1) ySol(1) zSol(1) xSol(2) ySol(2) zSol(2)]
    recalculate = 1; 
end

