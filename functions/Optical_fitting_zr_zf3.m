function [mub, mus, zr, zf, shift, Rsq] = Optical_fitting_zr_zf3(co, cross, datapath, mus_depth, sur_px)
% Fitting scattering and retardance of brain tissue
% co: volume data for co pol
% cross: volume data for cross pol
% s_seg: slice number
% Z_seg: tile number
% datapath: datapath for original data
% aip_threshold: intensity threshold to remove agarose, use same value as
% in OCT_recon.m
% mus_depth: number of z pixels for fitting mus
% bfg_depth: number of z pixels for fitting birefringence
% ds_factor: downsample factor
% zf: focus depth at center of FOV, i.e, x=middle, y=middle, starting from
% top of volume
% zf_tilt: zf difference along x axis of FOV

cd(datapath);
co1=imresize3(co,[size(co,1), size(co,2)/20, size(co,3)/20]);
cross1=imresize3(cross, [size(cross,1), size(cross,2)/20, size(cross,3)/20]);
ref=single(sqrt(co1.^2+cross1.^2));
sur=surprofile_agar(ref,'PSOCT',1); % find surface, downsample=20
Z_step=3;

% sensitivity roll-off correction
% w=2.2; % sensitivity roff-off constant, w=2.2 for 5x obj, w=2.22 for 10x obj
% I=rolloff_corr(I,w);    
vol=ref(:,:,:);
    %% Prefit with 2 parameters, zr = 70 um, zf = (zf-surface)*Z_step and adjust according to zf_tilt

N_alines = size(vol,2)*size(vol,3);

% Initial guess for the parameters
%initialGuess = [0.01, 0.01, 70, 70];  % Adjust the initial guess based on your model

%Reshape data into 2D vector
start_depth = sur+sur_px; %*** fix start depth to be a surface 2D
fit_depth=min(size(vol,1)-5,start_depth+mus_depth-1);%***
vol2 = vol(start_depth:fit_depth,:,:);
aline_len = fit_depth-start_depth;
for i=1:size(vol2,2)
    for j=1:size(vol2,3)
        z_1D = transpose(((0:aline_len(i,j))+sur_px)*Z_step); % *** make sure that zdata are custom for each aline
        z_3D(:,i,j) = z_1D;%repmat(z_1D,1,N_alines);
        x = 1:length(z_1D);
        w_3D(:,i,j) = x/length(x);
    end
end
%%
%     param = zeros(size(vol,2),size(vol,3),4); 

% Define the fitting function

% Generate sample data (replace this with your actual data)
ydata = double(vol2(:));  % Replace with your actual ydata matrix
%ydata = smooth(ydata);
zdata = z_3D(:);  % Replace with your actual zdata matrix
w = w_3D(:); %weights for fitting

% Initial guess for parameters (replace with your initial guess)
constant_initial_guess = [0.015, 0.1, 70, 0.003];
param4_initial_guess = 30 + 100*rand(N_alines, 1); % rand for param4

% Combine initial guesses
initial_guess = [constant_initial_guess, param4_initial_guess'];

% Set lower and upper bounds for parameters (if needed)
lower_bound = [0.00000001,  0.000001, 1,      0,      (1*ones(size(param4_initial_guess)))'];
upper_bound = [50,          50,    400,    1,    (400*ones(size(param4_initial_guess)))'];

% Perform simultaneous least squares fitting using lsqnonlin
options = optimoptions(@lsqnonlin, 'Display', 'iter');  % Adjust options as needed
% Create a function handle for the fitting function
fitting_function_handle = @(p) customFittingFunction(p, zdata, N_alines, aline_len(:));

% Perform fitting
[parameters, resnorm, residual] = lsqnonlin(@(p) w.*(ydata - fitting_function_handle(p)), initial_guess, lower_bound, upper_bound, options);

% Extract the constant parameters and the varying parameter
constant_parameters = parameters(1:4);
varying_parameter_matrix = parameters(5:N_alines+4);

% Reshape the varying parameter to match the original matrix size
varying_parameter_matrix = reshape(varying_parameter_matrix, size(vol,2), size(vol,3));
%%
% plot fitting results
% Extract parameters
constant_params = constant_parameters;
param4_matrix = parameters(5:N_alines+4);

% Initialize result
result = zeros(size(zdata));
aline_len = aline_len(:);
% Compute the result for each column in param4_matrix
for col = 1:length(param4_matrix)
        pxS = 1 + (aline_len(col)+1)*(col-1);
        pxE = pxS + (aline_len(col));
        out = double(constant_params(4)+sqrt(constant_params(1).* exp(-2 .* constant_params(2) .* zdata(pxS:pxE)) .* ...
            (1 ./ (1 + ((zdata(pxS:pxE) - (param4_matrix(col))) ./ constant_params(3)).^2))));
        result(pxS:pxE) = out;
        y = ydata(pxS:pxE);
        r2(col) = 1 - sum((y - out).^2)/sum((y - mean(y)).^2);
end
 figure(1); plot(result(1000:10500));hold on; plot(ydata(1000:10500));hold off; legend('fitting', 'data');
 zr=parameters(3);
 zf = varying_parameter_matrix;
 mean(mean(zf))
 mub = constant_params(1);
 mus = constant_params(2);
 shift = constant_params(4);
 Rsq = 1 - sum((ydata - result).^2)/sum((ydata - mean(ydata)).^2);
 r2_2D = reshape(r2, size(vol,2), size(vol,3));
 view3D(reshape(residual, size(vol2)));
 figure(2);imagesc(zf); colorbar; caxis([10 55]); title('Zf map');
 figure(3);imagesc(r2_2D);colorbar; caxis([0.8 1]); title('R2');
end



