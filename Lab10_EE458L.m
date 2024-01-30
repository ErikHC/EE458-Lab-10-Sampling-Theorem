% EE458L Lab 10
% Erik Chavarin

clear;
close all;

% Sampling freqs
Fs1 = 48e3;
Ts1 = Fs1^-1;
Fs2 = 8e3;
Ts2 = Fs2^-1;

% Signal freqs (Tasks 1 & 4)
fm1 = 2e3;
Tm1 = fm1^-1;
fm2 = 6e3;
Tm2 = fm2^-1;

% Time vector and resampled time vector
% Task 1:
t1 = 0:Ts1:Tm1*10-Ts1;
t_resamp1 = 0:Ts2:Tm1*10-Ts2;  % Needed to produce resampled signal elements
% Task 4:
t2 = 0:Ts1:Tm2*10-Ts1;
t_resamp2 = 0:Ts2:Tm2*10-Ts2; 
% Vectors must begin at 0 to create a fiducial marker

% Task 1: Signal generation
% fm = 2 kHz
x1 = sin(2 * pi * fm1 .* t1);
x_resamp1 = sin(2 * pi * fm1 .* t_resamp1);
% Task 4: fm = 6 kHz
x2 = sin(2 * pi * fm2 .* t2);
x_resamp2 = sin(2 * pi * fm2 .* t_resamp2);

% Task 2: Signal resampling
% Every [pad] zeros, an element from the resampled signal is inserted into 
% the vector so length(x) == length(y), same time vector may be used
pad = 6;  % T, 7T, 13T, ...
y1 = zeros(1,pad*length(x_resamp1));
y1(1:pad:end) = x_resamp1;
% Task 4:
y2 = zeros(1,pad*length(x_resamp2));
y2(1:pad:end) = x_resamp2;  % 80/6 \approx 13
% y2 = [y2, [0, 0]];  % Padding with two more zeros

% Task 3: LPF applied to resampled signal
%firpm LPF specs: 4 kHz pb + 400 Hz tb, Fs = 48 kHz
pb = 4e3;  % passband width
tb = 4e2;  % transition bandwidth
fI = [0, pb, pb+tb, Fs1/2] / (Fs1 / 2);  % Normalized by Nyquist freq
aI = [1, 1, 0, 0];  % LPF with full suppression on output
N = 31;  % Order 31, firpm returns N+1 samples
g = firpm(N, fI, aI);  % Impulse response
G = fftshift(fft(g));  % Transfer function centered about DC
fshift = (-N/2:N/2)*(Fs1/N);  % Symmetric spectrum axis


z1 = conv(g, y1);
z2 = conv(g, y2);
t_conv1 = 0:Ts1:(Tm1*10+Ts1*N)-Ts1;
t_conv2 = 0:Ts1:(Tm2*10+Ts1*N)-3*Ts1;  % Shave off a bit (because 80/6)
% t_conv2 = 0:Ts1:(Tm2*10+Ts1*N)-Ts1;  % N = 80 (no shaving)

% Plotting (Signals)
sgtitle("Lab 10: Nyquist--Shannon Sampling Theorem", ...
    "FontSize", 14, "FontName", "Serif", "Interpreter", "latex")
% x1
subplot(2,3,1)
plot(t1*1e3, x1)
title("Task 1: Unsampled Signal $x_1(t) = \sin(2 \pi t \cdot 2000)$, " + ...
    "$F_{s1} = 48$kHz", "FontSize", 12, "FontName", "Serif", "Interpreter", ...
    "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$x_1(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
ylim([-1.5,1.5]);

subplot(2,3,2)
stem(t1*1e3, y1)
title("Task 2: Resampled Signal, $y_1(t)$, $F_{s2} = 8$kHz ", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$y_1(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
ylim([-1.5,1.5]);

subplot(2,3,3)
plot(t_conv1*1e3, z1)
title("Task 3: $z_1(t) = g(t) * y_1(t)$, PB = 4kHz, $N_{LPF}$ = 31", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$z_1(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")

% x2
subplot(2,3,4)
plot(t2*1e3, x2)
title("Task 1: Unsampled Signal $x_2(t) = \sin(2 \pi t \cdot 6000)$, " + ...
"$F_{s1} = 48$kHz", "FontSize", 12, "FontName", "Serif", "Interpreter", ...
    "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$x_2(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
ylim([-1.5,1.5]);

subplot(2,3,5)
stem(t2(1:1:length(t2)-2)*1e3, y2)
% stem(t2*1e3, y2)
title("Task 2: Resampled Signal, $y_2(t)$, $F_{s2} = 8$kHz ", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$y_2(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
ylim([-1.5,1.5]);

subplot(2,3,6)
plot(t_conv2*1e3, z2)
title("Task 3: $z_2(t) = g(t) * y_2(t)$, PB = 4kHz, $N_{LPF} = 31$ ", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xlabel("t (ms)", "FontSize", 12, "FontName", "Serif")
ylabel("$z_2(t)$", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")

% Plotting (LPF)
plot(fshift/1000, 20*log10(abs(G)))
title("Task 3: Power Spectrum of $4$kHz Low-pass Reconstruction Filter", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xlabel("$f$ (kHz)", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
ylabel("$|G(f)|^2$ (dB)", ...
    "FontSize", 12, "FontName", "Serif", "Interpreter", "latex")
xticks([-24 -4 0 4 24])  % Plots frequency points of interest
xline([-Fs1/2e3, Fs1/2e3], "-.", ["-Fs/2", "Fs/2"], ...
    LabelHorizontalAlignment="center", LabelVerticalAlignment="middle")

% Label passband and transition band
xline(pb/1e3, "-.", "4 kHz", ...
    LabelHorizontalAlignment="left", LabelVerticalAlignment="bottom")
xline( (pb+tb)/1e3, "-.", "4.4 kHz", ...
    LabelHorizontalAlignment="right", LabelVerticalAlignment="bottom")
ylim([-40, 10])
