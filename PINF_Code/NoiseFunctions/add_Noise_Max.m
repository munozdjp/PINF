function [xNoisy, Noise_normalized]=add_Noise_Max(x,noise_scale)
%if I scale with no normalization i need to uncoment the normalization
%funciton after the integration in the main code

%Noise for the scaling with the maximum 
%     ss = size(x);
%     Noise_scaled=noise_scale*randn(ss);% Noise N(0,scale)
%     maximum = max(abs(x));% maximum on each stat variable
%     Noise_normalized = Noise_scaled.*maximum; %Normalized with respect maximum
%     xNoisy= x + Noise_normalized;

%Noise for the scaling respect the variable x, Noise = x*noise
%     ss = size(x);
%     Noise_scaled=noise_scale*randn(ss);% Noise N(0,scale)
%     Noise_normalized = Noise_scaled.*x; %Normalized with respect maximum
%     xNoisy = x + Noise_normalized;

%If I want to use the Noise-Normalized I need to do Normalize the variable 
%"Command of normalizing is after the integration"
%Noise when I normalized my variables x
    ss = size(x);
    Noise_scaled =noise_scale*randn(ss);% Noise N(0,scale)
    Noise_normalized = Noise_scaled;
    xNoisy= x + Noise_scaled;
end