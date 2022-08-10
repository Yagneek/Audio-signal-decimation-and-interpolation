clear all;
close all;
clc;

[d,r] = audioread('msmn1.wav');
figure('Name', 'Spectrum of the input sound file')
specgram(d,1024,r)
soundsc(d,r)
pause(size(d)/r);

for i = 1:3
    M=2^i;
    if i==1
        Q = 29;
    elseif i==2
        Q = 28;
    else
        Q = 27;
    end
    d_decimated = fxd_decimator(d,M,r/(2*M),r,1001,1,Q);
    figure('Name', ['Spectrum after decimation by a factor of ' int2str(2^i)])
    specgram(d_decimated,1024,r/M)
    soundsc(d_decimated,r/M)
    audiowrite(strcat('decimated_',int2str(2^i),'.wav'),d_decimated,floor(r/M));
    pause(size(d_decimated)/(r/M));
    
    d_interpolated = fxd_interpolator(d_decimated,M,r/(2*M),r,1001,M,Q);
    figure('Name', ['Spectrum after decimation and interpolation by a factor of ' int2str(2^i)])
    specgram(d_interpolated,1024,r)
    soundsc(d_interpolated,r)
    audiowrite(strcat('interpolated_',int2str(2^i),'.wav'),d_interpolated,r);
    pause(size(d_interpolated)/r);
end

function w = hamming_window(N)  %generates a hamming window with N samples
    w = zeros(1,N);
    for i=1:N
        n = i-1;    %since the array elements are indexed from 1
        w(i) = 0.54 - (0.46*cos(2*pi*(n/(N-1))));
    end
end

%outputs the h[n] of a low pass filter with cutoff frequency w_c(obtained from fc and fs) and N samples
function h = lpf(fc,fs,N,gain)
    h_d = zeros(1,N);
    w_c = 2*pi*(fc/fs);
    for i=1:N
        n = i-((N+1)/2);    %since the array elements are indexed from 1
        if n ~= 0
            h_d(i) = (sin(w_c*n))/(pi*n);
        else
            h_d(i) = w_c/pi;
        end
    end
    W_H = hamming_window(N);
    h = h_d.*W_H;
    h = h*gain;
end

function h = fxd_lpf(fc,fs,N,gain,Q)
    h_d = zeros(1,N);
    w_c = 2*pi*(fc/fs);
    for i=1:N
        n = i-((N+1)/2);    %since the array elements are indexed from 1
        if n ~= 0
            h_d(i) = (sin(w_c*n))/(pi*n);
        else
            h_d(i) = w_c/pi;
        end
    end
    W_H = hamming_window(N);
    fxd_h_d = toFixed(h_d,Q);
    fxd_W_H = toFixed(W_H,Q);
    fxd_h = fxdMul(fxd_h_d,Q,fxd_W_H,Q,Q);
    fxd_gain = toFixed(gain,5);
    fxd_h = fxdMul(fxd_h,Q,fxd_gain,5,Q);
    h = toFloat(fxd_h,Q);
end

%Downsamples the given input signal x[n] by a factor of M
function x_down = downsampler(x,M)
    lx = length(x);
    x_down = zeros(floor((lx-1)/M)+1,1);
    for i=1:length(x_down)
        x_down(i) = x((M*(i-1))+1); %Since the array elements are indexed from 1
    end
end

%Decimates the given input signal x[n] by a factor of M with 
%anti-aliasing filter specifications fc,fs,N=lh and given gain
function xd = decimator(x,M,fc,fs,lh,gain)
    anti_aliasing = lpf(fc,fs,lh,gain);
    xf = conv(x,anti_aliasing); %Passing the input signal x[n] through the anti-aliasing filter
    lx = length(x);
    xf = xf((lh-1)/2+1:(lh-1)/2+lx);    %Taking only the middle lx samples
    xd = downsampler(xf,M);
end

function xd = fxd_decimator(x,M,fc,fs,lh,gain,Q)
    anti_aliasing = fxd_lpf(fc,fs,lh,gain,Q);
    fxd_x = toFixed(x,Q);
    fxd_anti_aliasing = toFixed(anti_aliasing,Q);
    fxd_xf = fxd_conv(fxd_x,fxd_anti_aliasing,Q); %Passing the input signal x[n] through the anti-aliasing filter
    lx = length(fxd_x);
    fxd_xf = fxd_xf((lh-1)/2+1:(lh-1)/2+lx);    %Taking only the middle lx samples
    xf = toFloat(fxd_xf,Q);
    xd = downsampler(xf,M);
end

%Upsamples the given input signal x[n] by a factor of L
function xu = upsampler(x,L)
    lx = length(x);
    xu = zeros(L*lx,1);
    for i=1:length(xu)
        if mod(i-1,L)==0
            xu(i) = x((i-1)/L+1);   %Since the array elements are indexed from 1
        else
            xu(i) = 0;
        end
    end
end

%Interpolates the given input signal x[n] by a factor of L with 
%anti-imaging filter specifications fc,fs,N=lh and given gain
function xi = interpolator(x,L,fc,fs,lh,gain)
    xu = upsampler(x,L);
    anti_imaging = lpf(fc,fs,lh,gain);
    xi = conv(xu,anti_imaging); %Passing the upsampled signal xu[n] through the anti-imaging filter
    lxu = length(xu);
    xi = xi((lh-1)/2+1:(lh-1)/2+lxu);  %Taking only the middle lxu samples
end

function xi = fxd_interpolator(x,L,fc,fs,lh,gain,Q)
    xu = upsampler(x,L);
    anti_imaging = fxd_lpf(fc,fs,lh,gain,Q);
    fxd_xu = toFixed(xu,Q);
    fxd_anti_imaging = toFixed(anti_imaging,Q);
    fxd_xi = fxd_conv(fxd_xu,fxd_anti_imaging,Q); %Passing the upsampled signal xu[n] through the anti-imaging filter
    lxu = length(xu);
    fxd_xi = fxd_xi((lh-1)/2+1:(lh-1)/2+lxu);  %Taking only the middle lxu samples
    xi = toFloat(fxd_xi,Q);
end

function res = fxdMul(num1,q1,num2,q2,resq)
    res = double(num1).*double(num2);
    res = int32(res/2^(q1+q2-resq));
end

function output = toFixed(A,Q)
    output = A.*2^Q;
    output = int32(output);
end

function output = toFloat(A,Q)
    A = double(A);
    output = A./2^Q;
end

function res = fxdAdd(A,B)
    res = A+B;
end

function Y = fxd_conv(X,H,Q)
    len_H = length(H);
    len_X = length(X);
    len_Y = len_X+len_H-1;
    Y = int32(zeros(len_Y,1));
    for n=1:len_Y
        for k=1:n
            if n-k+1>len_H||k>len_X
                continue;
            end
            Y(n) = fxdAdd(Y(n),fxdMul(X(k),Q,H(n-k+1),Q,Q));
        end
    end
end