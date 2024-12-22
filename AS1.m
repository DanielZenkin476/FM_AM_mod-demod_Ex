
close all

[v_m, fs] = audioread('in-the-air.wav') ; % Load the WAV file into v_m(signal) and fs(sample rate)
v_m = v_m(:,1); %using one channel only
%sound(v_m,fs);% play wav file to check if loaded correctly
%defining parametars
Ts = 1/fs;
N = length(v_m);
t = 0:Ts:(N-1)*Ts;
f = linspace(-fs/2,fs/2,N);
%compute fft on signal
V_m = fftshift(fft(v_m)) / sqrt(N);

%1.1.1 -Plot on the same figure (using “subplot”) the data signal and its Fourier transform 
figure;
subplot(2,1,1);% 1st subplot for v_m(t)
plot(t,v_m);
xlabel('Time (s)')
ylabel('v_m(t)')
title("v_m(t)")
grid on

subplot(2,1,2);% 2st subplot for V_m(F)
plot(f,abs(V_m));
xlabel('Freq. (Hz)')
ylabel('|V_m(f)|')
title("V_m(f)")
grid on


%1.2

fc = 15*10^3; %carrier frequency
K_AM = 0.02;% modulation index

%1.2.1Create the signal v_AM

v_AM = ammod(v_m, fc, fs, 0, K_AM);% AM mod with carrier amp. K_AM

%1.2.2: Compute the Fourier transform of v_AM(t) and plot it 

V_AM = fftshift(fft(v_AM)) / sqrt(N);  % Fourier Transform normalized


figure;
plot(f, abs(V_AM));  
xlabel('Frequency (Hz)');
ylabel('|V_AM(f)|');
title('|V_AM(f)|');
grid on 


%1.2.3
%sound(v_AM,fs); % we dont hear anything- signal modulated on 15 khz carier which humans cant hear 
%1.3.1. Generate in MATLAB a real white Gaussian noise sequence

z = 0.02 * randn(1,N);% sqrt(n_0/2) = 0.02 
z = reshape(z,size(v_AM));% size of z was not matching so reshaped it so an error wouldnt jump- array was 89.13 Gb - matlab said it couldnt handle it)

%1.3.2 Generate the received signal at the channel output. Compute and plot in MATLAB the Fourier transform

%generate xr(t)
x_r = v_AM + z;
%generate Xr(f)
X_r = fftshift(fft(x_r)) / sqrt(N);

%Plot the time-domain signal x_r(t)
figure;
subplot(2, 1, 1);
plot(t, x_r);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal x_r(t)');
grid on

%Plot the frequency-domain signal X_r(f)
subplot(2, 1, 2);
plot(f, abs(X_r));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_r(f)|');
title('Fourier Transform |X_r(f)|');
grid on 

%1.4
fm = 5*10^3;
x_L = bandpass(x_r, [fc-fm fc+fm], fs);  

% Compute the Fourier Transform of x_L(t)
X_L = fftshift(fft(x_L)) / sqrt(N);  % Fourier Transform normalized

% Plot the filtered signal x_L(t) and its Fourier Transform
figure;
subplot(2, 1, 1);
plot(t, x_L);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal x_L(t)');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_L));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_L(f)|');
title('Fourier Transform |X_L(f)|');
grid on

%1.4.2

x_d = amdemod(x_L,fc,fs,0,K_AM);
x_d = lowpass(x_d,fm,fs);
X_d = fftshift(fft(x_d)) / sqrt(N);

figure;
subplot(2, 1, 1);
plot(t, x_d);
hold on
plot(t,v_m)
xlabel('Time (s)');
ylabel('Amplitude');
legend('x_d(t)','v_m(t)')
title('x_d(t) and V_m(t)');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_d));  % Magnitude of Fourier Transform
hold on
plot(f,(abs(V_m)))
xlabel('Frequency (Hz)');
ylabel('|signal(f)|');
legend('|X_d(f)|','|V_m(f)|')
title('|X_d(f)| and |V_m(f)|');
grid on
%1.4.3
%sound(x_d,fs)

%1.4.4

disp('Correlation x_d and v_m, (sqrt(N_0/2)=0.02)  :')
disp(xcorr(x_d, v_m, 0, 'coeff'))

%1.5
% with sqrt(N_0/2) = 0.1

%1.3.1. Generate in MATLAB a real white Gaussian noise sequence

z_2 = 0.1* randn(1,N);% sqrt(n_0/2) = 0.01
z_2 = reshape(z_2,size(v_AM));% size of z was not matching so reshaped it so an error wouldnt jump


%1.3.2 Generate the received signal at the channel output. Compute and plot in MATLAB the Fourier transform

%generate xr(t)
x_r_2 = v_AM + z_2;
%generate Xr(f)
X_r_2 = fftshift(fft(x_r_2)) / sqrt(N);

%Plot the time-domain signal x_r(t)
figure;
subplot(2, 1, 1);
plot(t, x_r_2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal x_r(t),sqrt(N_0/2) = 0.1');
grid on

%Plot the frequency-domain signal X_r(f)
subplot(2, 1, 2);
plot(f, abs(X_r_2));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_r(f)|');
title('Fourier Transform |X_r(f)|,sqrt(N_0/2) = 0.1');
grid on 

%The plot shows the frequency-domain representation of the received signal.
% It includes the original signal components at the carrier and sidebands,
% along with noise spread across all frequencies.

%1.4
x_L_2 = bandpass(x_r_2, [fc-fm fc+fm], fs); 

% Compute the Fourier Transform of x_L(t)
X_L_2 = fftshift(fft(x_L_2)) / sqrt(N);  % Fourier Transform normalized

% Plot the filtered signal x_L(t) and its Fourier Transform
figure;
subplot(2, 1, 1);
plot(t, x_L_2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal x_L(t),sqrt(N_0/2) = 0.1');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_L_2));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_L(f)|');
title('Fourier Transform |X_L(f)|,sqrt(N_0/2) = 0.1');
grid on

% Explanation:
% - The time-domain plot shows the filtered signal around the carrier frequency.
% - The frequency-domain plot highlights that the noise outside the bandpass range has been attenuated, leaving the carrier and sidebands.

%1.4.2

x_d_2 = amdemod(x_L_2,fc,fs,0,K_AM);%demod x_l
x_d_2 = lowpass(x_d_2,fm,fs);% pass signal thought a lowpass filter
X_d_2 = fftshift(fft(x_d_2)/sqrt(N));% calc. foriuer transform

figure;
subplot(2, 1, 1);
plot(t, x_d_2);
hold on
plot(t,v_m)
xlabel('Time (s)');
ylabel('Amplitude');
legend('x_d(t)','v_m(t)')
title('x_d(t) and V_m(t),sqrt(N_0/2) = 0.1');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_d_2));  % Magnitude of Fourier Transform
hold on
plot(f,(abs(V_m)))
xlabel('Frequency (Hz)');
ylabel('|signal(f)|');
legend('|X_d(f)|','|V_m(f)|')
title('|X_d(f)| and |V_m(f)|,sqrt(N_0/2) = 0.1');
grid on
%1.4.3
%sound(x_d_2,fs)

%1.4.4

disp('Correlation x_d and v_m(sqrt(N_0/2) = 0.1) :')
disp(xcorr(x_d_2, v_m, 0, 'coeff'))


%1.5.2 -Fm modulation 
d_fd = 10*10^3;
K_FM = d_fd/max(v_m);
v_FM = fmmod(v_m,fc,fs,d_fd,0);
V_FM = fftshift(fft(v_m)/sqrt(N));

%sound(v_FM,fs);

%1.3.1. Generate in MATLAB a real white Gaussian noise sequence


%1.3.2 Generate the received signal at the channel output. Compute and plot in MATLAB the Fourier transform

%generate xr(t)
x_r_FM = v_FM + z;% using old z 
%generate Xr(f)
X_r_FM = fftshift(fft(x_r_FM)) / sqrt(N);

%Plot the time-domain signal x_r(t)
figure;
subplot(2, 1, 1);
plot(t, x_r_FM);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal x_r(t), FM sqrt(n_0/2) = 0.02 ');
grid on

%Plot the frequency-domain signal X_r(f)
subplot(2, 1, 2);
plot(f, abs(X_r_FM));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_r(t)');
title('Fourier Transform |X_r(f)|,FM sqrt(n_0/2) = 0.02');
grid on 

%The plot shows the frequency-domain representation of the received signal.
% It includes the original signal components at the carrier and sidebands,
% along with noise spread across all frequencies.

%1.4
x_L_FM = bandpass(x_r_FM, [fc-fm fc+fm], fs); 

% Compute the Fourier Transform of x_L(t)
X_L_FM = fftshift(fft(x_L_FM)) / sqrt(N);  % Fourier Transform normalized

% Plot the filtered signal x_L(t) and its Fourier Transform
figure;
subplot(2, 1, 1);
plot(t, x_L_FM);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal x_L(t),FM sqrt(n_0/2) = 0.02');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_L_FM));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_L(f)|');
title('Fourier Transform |X_L(f)|,FM sqrt(n_0/2) = 0.02');
grid on

%1.4.2

x_d_FM = fmdemod(x_L_FM,fc,fs,d_fd,0);
x_d_FM = lowpass(x_d_FM,fm,fs);
X_d_FM = fftshift(fft(x_d_FM)) /sqrt(N);

figure;
subplot(2, 1, 1);
plot(t, x_d_FM);
hold on
plot(t,v_m)
xlabel('Time (s)');
ylabel('Amplitude');
legend('x_d(t)','v_m(t)')
title('x_d(t) and V_m(t)');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_d_FM));  % Magnitude of Fourier Transform
hold on
plot(f,(abs(V_m)))
xlabel('Frequency (Hz)');
ylabel('|signal(f)|');
legend('|X_d(f)|','|V_m(f)|')
title('|X_d(f)| and |V_m(f)|');
grid on
%1.4.3
%sound(x_d_FM,fs)

%1.4.4

disp('Correlation x_d and v_m :,FM sqrt(n_0/2) = 0.02')
disp(xcorr(x_d_FM, v_m, 0, 'coeff'))

%1.5
% with sqrt(N_0/2) = 0.1

%1.3.1. Generate in MATLAB a real white Gaussian noise sequence


%1.3.2 Generate the received signal at the channel output. Compute and plot in MATLAB the Fourier transform

%generate xr(t)
x_r_FM_2 = v_FM + z_2;
%generate Xr(f)
X_r_FM_2 = fftshift(fft(x_r_FM_2)) / sqrt(N);

%Plot the time-domain signal x_r(t)
figure;
subplot(2, 1, 1);
plot(t, x_r_FM_2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal x_r(t), FM sqrt(N_0/2) = 0.1');
grid on

%Plot the frequency-domain signal X_r(f)
subplot(2, 1, 2);
plot(f, abs(X_r_FM_2));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_r(t)|');
title('Fourier Transform |X_r(f)|,FM sqrt(N_0/2) = 0.1');
grid on 

%1.4

x_L_FM_2 = bandpass(x_r_FM_2, [fc-fm fc+fm], fs); 

% Compute the Fourier Transform of x_L(t)
X_L_FM_2 = fftshift(fft(x_L_FM_2)) / sqrt(N);  % Fourier Transform normalized

% Plot the filtered signal x_L(t) and its Fourier Transform
figure;
subplot(2, 1, 1);
plot(t, x_L_FM_2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass Filtered Signal x_L(t),FM sqrt(N_0/2) = 0.1');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_L_FM_2));  % Magnitude of Fourier Transform
xlabel('Frequency (Hz)');
ylabel('|X_L(f)|');
title('Fourier Transform |X_L(f)|,FM sqrt(N_0/2) = 0.01');
grid on

%1.4.2

x_d_FM_2 = fmdemod(x_L_FM_2,fc,fs,d_fd,0);
x_d_FM_2 = lowpass(x_d_FM_2,fm,fs);
X_d_FM_2 = fftshift(fft(x_d_FM_2)/sqrt(N));

figure;
subplot(2, 1, 1);
plot(t, x_d_FM_2);
hold on
plot(t,v_m)
xlabel('Time (s)');
ylabel('Amplitude');
legend('x_d(t)','v_m(t)')
title('x_d(t) and V_m(t)');
grid on 

subplot(2, 1, 2);
plot(f, abs(X_d_FM_2));  % Magnitude of Fourier Transform
hold on
plot(f,(abs(V_m)))
xlabel('Frequency (Hz)');
ylabel('|signal(f)|');
legend('|X_d(f)|','|V_m(f)|')
title('|X_d(f)| and |V_m(f)|');
grid on
%1.4.3
sound(x_d_FM_2,fs)

%1.4.4

disp('Correlation x_d and v_m(FM sqrt(N_0/2) = 0.01) :')
disp(xcorr(x_d_FM_2, v_m, 0, 'coeff'))








