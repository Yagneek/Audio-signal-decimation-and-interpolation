clear all;
close all;
clc;
%The given input signal
lx = 96;
f0 = 100;
f1 = 200;
f2 = 300;
fs = 2400;
x = zeros(1,lx);
for i=1:lx
    n = i-1;
    x(i) = sin(2*pi*f0*n/fs)+0.5*sin(2*pi*f1*n/fs)+0.6*sin(2*pi*f2*n/fs);
end
%fvtool(x,1) %Plotting the x[n]
%figure(1);
%stem(x)
%Decimation and interpolation by factor 2 with the given filter
%specifications : M=L=2;fc=600Hz;fs=2400Hz;N=101
%Anti-aliasing gain=1;Anti-imaging gain=L
for i=1:3
    M=2^i;
    if i==1
        Q = 29;
    elseif i==2
        Q = 28;
    else
        Q = 27;
    end
    xd1 = decimator(x,M,600,2400,101,1);
    fxd_xd1 = fxd_decimator(x,M,600,2400,101,1,Q);
    % fvtool(xd1,1)   %Plotting xd[n]
    y1 = interpolator(xd1,M,600,2400,101,M);
    fxd_y1 = fxd_interpolator(fxd_xd1,M,600,2400,101,M,Q);
    % fvtool(y1,1)    %Plotting y[n]
    e1 = fxd_y1-y1;  %Computing the error vector
    figure();
    stem(fxd_y1)
    figure();
    stem(e1)    %Plotting the error vector
    %Calculating the average error
    msq_error1 = 0;
    for i=1:lx
        msq_error1 = msq_error1 + (e1(i)^2);
    end
    msq_error1 = msq_error1/lx;
    msq_error1 = sqrt(msq_error1)
end

% a = int32([1,2,5,7,3,4,8]);
% b = int32([4,5,7,3]);
% o = fxd_conv(a,b,4)


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
    x_down = zeros(1,floor((lx-1)/M)+1);
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
    xf = xf((lh-1)/2+1:(lh-1)/2+lx);   %Taking only the middle lx samples
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
    %max_xf = max(xf)
    xd = downsampler(xf,M);
end

%Upsamples the given input signal x[n] by a factor of L
function xu = upsampler(x,L)
    lx = length(x);
    xu = zeros(1,L*lx);
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
    %max_lpf = max(anti_imaging)
    fxd_xu = toFixed(xu,Q);
    fxd_anti_imaging = toFixed(anti_imaging,Q);
    fxd_xi = fxd_conv(fxd_xu,fxd_anti_imaging,Q); %Passing the upsampled signal xu[n] through the anti-imaging filter
    lxu = length(xu);
    fxd_xi = fxd_xi((lh-1)/2+1:(lh-1)/2+lxu);  %Taking only the middle lxu samples
    xi = toFloat(fxd_xi,Q);
    %max_xi = max(xi)
end

function res = fxdMul(num1,q1,num2,q2,resq)
    num1 = cast(num1,"double");
    num2 = cast(num2,"double");
    res = (num1.*num2)/bitshift(1,q1+q2-resq);
    res = cast(res,"int32");
end

function output = toFixed(A,Q)
    output = A.*bitshift(1,Q);
    output = cast(output,"int32");
end

function output = toFloat(A,Q)
    A = cast(A,"double");
    output = A./bitshift(1,Q);
end

function res = fxdAdd(A,B)
    res = A+B;
end

function Y = fxd_conv(X,H,Q)
    Y = int32(zeros(1,length(X)+length(H)-1));
    for n=1:length(Y)
        for k=1:n
            if n-k+1>length(H)||k>length(X)
                continue;
            end
            Y(n) = fxdAdd(Y(n),fxdMul(X(k),Q,H(n-k+1),Q,Q));
        end
    end
end