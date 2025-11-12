clear all;close all;clc;
j=1;
%for fit_depth = 180:-1:50
fit_depth = 130;
N_iter = 10^3;

opts = optimset('Display','off','TolFun',1e-10);

Z_step=3; % in mm

C = 0.0285;
zr_arr = 43;%(58-43).*rand(N_iter,1)+43;
max_noise = 0.1;
mus_arr = 0.040.*rand(N_iter,1)+0.0005; 
mub_arr = 0.65;%rand(N_iter,1)+.5;
zf_arr = 130;%(150-50).*rand(N_iter,1) + 50;


for i = 1:N_iter
    zdata = (1:fit_depth).*Z_step;
    mus = mus_arr(i); 
    %zf_noise = 50*rand(1,1);
    %zr_noise = 10*rand(1,1);
    zf = zf_arr;%(i);
    mub = mub_arr;%(i); 
    zr = zr_arr;%(i);

    %generate noise
    gamma_noise = randg(1,[1, fit_depth]);
    temp = movmean(gamma_noise,5); % slow down the noise to match real data
    noise = max_noise*temp;
    
    %create signal
    aline = double(C + sqrt(mub.*exp(-2.*mus.*(zdata)).*(1./(1+((zdata-zf)./zr).^2)))); 
    %add weighted noise to the signal
    aline_weight = aline./max(aline);
    noise_w = noise.*(aline_weight.^2);
    noise_arr(i,:) = noise_w;
    ref = aline + noise_w;
    ydata = ref;
    % define fitting function
    %zr = zr_arr(i)+zr_noise;
    fun_pix = @(p,zdata)double(C + sqrt(p(1).*exp(-2.*p(2).*(zdata)).*(1./(1+((zdata-zf)./zr).^2))));  
    % fitting
    lb = [10^(-10) 10^(-10) ];
    ub=[50 50 ];

    [param(i,1:2), rn] = lsqcurvefit(fun_pix,[0.01 0.01 ],zdata,ydata,lb,ub,opts);
    R_2(i) = 1 - sum((ydata - fun_pix(param(i,:),zdata)).^2)/sum((ydata - mean(ydata)).^2);
% 
    %figure(1); scatter(zdata,ydata);hold on;y=fun_pix(param(i,:),zdata);plot(zdata,y); %ylim([0 0.1]);        
end
figure(2);plot(500.*mus_arr, 500.*param(:,2), 'o');
est = param(:,2); real = mus_arr;
errors = 100.*(est-real)./real;
xlabel('theoretical mus');
ylabel('fitted, mus');
[Xsorted,I] = sort(real);
Ysorted = errors(I);
h = findobj(gca, 'Type', 'line');
set(h, 'LineWidth', 2);
% Make fonts bigger
set(gca, 'FontSize', 16);

figure(3); scatter(500.*Xsorted, abs(Ysorted)); ylabel('error, %'); 
xlabel('\mu_s, mm^{-1}');
err_arr(j) = std(errors); 
j=j+1;
clear noise_arr
%end
% Make all lines thicker
h = findobj(gca, 'Type', 'line');
set(h, 'LineWidth', 2);
% Make fonts bigger
xlim([0 20]);

set(gca, 'YScale', 'log');  % Make y-axis logarithmic
ylim([1e-3, 1e3]);
box off;
set(gca, 'FontSize', 20);

figure(4); hist(Ysorted, -30:10:100); xlabel('error, %'); 
%end
% Make all lines thicker
h = findobj(gca, 'Type', 'line');
set(h, 'LineWidth', 2);
% Make fonts bigger
set(gca, 'FontSize', 20);
xlim([-30 100])
box off;

error_std = std(errors)
error_mean = mean(errors)
error_median = median(errors)
error_range = [min(errors), max(errors)]