% clear all;close all;clc;
id='1'; %comment it when run fitting in parallel
% fitting scattering and birefringence of brain tissue
% run this after OCT_recon.m is finished
% clear all;close all;clc;
% read variable across samples parameters from .txt file
info = readlines('step0_info_d1.txt');
folder = info(2);      % OCT file path                          
P2path = info(5);   % 2P file path                          
datapath=strcat(folder,'dist_corrected/'); 
nslice=str2num(info(8)); % define total number of slices                                                                       
stitch=str2num(info(11)); % 1 means using OCT data to generate stitching coordinates, 
% 0 means using 2P stitching coordinates.                                                                        
ds_factor=str2num(info(14));     % downsampling factor, 4 means 4x4 pixel downsample, which is 12x12um pixel size                

% add subfunctions for the script. Change directory if not running on BU SCC
addpath(strcat(pwd, '\functions'));

cd(datapath);
create_dir(nslice, folder); 
sys = 'PSOCT';
% specify mosaic parameters, you can get it from Imagej stitching
xy=0;     % xy is the Y displacement of two adjacent tile align in the X direction, default to 0
yx=0;      % xx is the X displacement of two adjacent tile align in the Y direction, default to 0
numX=str2num(info(17));    % #tiles in X direction                                                                             
numY=str2num(info(20));    % #tiles in Y direction                                                                              

res=str2num(info(23));
if res==1 % if tile size 1x1 mm use following parameters
    xx=543;%900 866;    % xx is the X displacement of two adjacent tile align in the X direction
    yy=543;%900 866;    % yy is the Y displacement of two adjacent tile align in the Y direction
    pxlsize=[700 700];%[1000 1000];
    Xoverlap=0.15;%0.05;   % overlap in X direction
    Yoverlap=0.15;%0.05;   % overlap in Y direction

elseif res==3 % if tile size 3.3x3.3 mm use following parameters
    xx=866;
    yy=866;
    pxlsize=[1000 1000];
    Xoverlap=0.05;
    Yoverlap=0.05;
end

disp=[xx xy yy yx];
mosaic=[numX numY Xoverlap Yoverlap];
pattern = 'bidirectional';  % mosaic pattern, could be bidirectional or unidirectional
ntile=numX*numY;                                                                                                 
njobs=1;
section=ceil(ntile/njobs);
% the $SGE-TASK-ID environment variable read in is CHARACTER, need to transfer to number
id=str2num(id); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
istart=1;%!!!!
istop=section;
%% find by fitting averaged volume zf-map and zr-const
tic
for islice = 1%:floor(nslice/4):nslice % make a step size so we have 5 slices equally spaced across volume for averaging agar
    cd(datapath);
    filename0=dir(strcat('co-',num2str(islice),'-*.dat')); 
    name=strsplit(filename0(1).name,'.');  
    name_dat=strsplit(name{1},'-');
    slice_index=islice;
    %Xrpt and Yrpt are x and y scan repetition, default = 1
    Zsize = str2num(name_dat{4}); Xrpt = 1; Xsize=str2num(name_dat{5}); Yrpt = 1; Ysize = str2num(name_dat{6});
    mask_all=uint8(zeros(Xsize,Ysize));
    co_agar = zeros(Zsize,Xsize,Ysize); 
    cross_agar = zeros(Zsize,Xsize,Ysize); 
    agar = zeros(Xsize,Ysize);
    i=0;
   %-----------------------------------------------------------------------
    % !!!hand pick which tiles to use to find zf and zr (scattering agar) !!!!
    %----------------------------------------------------------------------
    for iFile = [1:690]
        name=strsplit(filename0(1).name,'.');  
        name_dat=strsplit(name{1},'-');
        slice_index=islice;
        %Xrpt and Yrpt are x and y scan repetition, default = 1
        Zsize = str2num(name_dat{4}); Xrpt = 1; Xsize=str2num(name_dat{5}); Yrpt = 1; Ysize = str2num(name_dat{6});
        dim1=[Zsize Xrpt Xsize Yrpt Ysize];     % tile size for reflectivity 
        name1=strcat(datapath,'co-',num2str(islice),'-',num2str(iFile),'-',num2str(Zsize),'-',num2str(Xsize),'-',num2str(Ysize),'.dat'); % gen file name for reflectivity
        if isfile(name1)
            %load reflectivity data
            co = ReadDat_int16(name1, dim1)./65535*4; 
            name1=strcat(datapath,'cross-',num2str(islice),'-',num2str(iFile),'-',num2str(Zsize),'-',num2str(Xsize),'-',num2str(Ysize),'.dat'); % gen file name for reflectivity
            cross = ReadDat_int16(name1, dim1)./65535*4; 
        else
            co=zeros(Zsize,Xsize,Ysize);
            cross=zeros(Zsize,Xsize,Ysize);
        end
    
        message=strcat('Tile No. ',string(iFile),' is read.', datestr(now,'DD:HH:MM'),'\n');
        fprintf(message);
    

        %save all the pixels
        ref=single(sqrt(co.^2+cross.^2));
        aip=squeeze(mean(ref(1:110,:,:),1));
        % 
        %figure();
        %imagesc(aip);
        %title(num2str(iFile));
        %view3D(ref);

        %save data for further averaging
          i = i+1;

         co_agar = co_agar+co;
         cross_agar = cross_agar+cross;
    
    end
    co = fillmissing(co_agar, 'nearest', 2);
    cross = fillmissing(cross_agar, 'nearest', 2);
    co = co./i; cross = cross./i; %average tiles
    cd(strcat(folder,'fitting/'));
    save(strcat('ave_tile',num2str(islice),'.mat'), "co", "cross");
    
end
toc

%% find by fitting average tile zf-map and zr-const
[mub, mus,zr, zf2, shift] = Optical_fitting_zr_zf3(co, cross, folder, 130, 5);
cd(strcat(folder,'fitting/'));
[X, Y] = meshgrid(1:35, 1:35);
[fitresult, gof] = createPolyFitOrder3(X, Y, zf2);
zf_l = fitresult(X,Y);
save('zf_zr_sim_fit.mat','zr', 'zf2', 'zf_l', 'mub', 'mus', 'shift');

%% plot  fitting results
co1=imresize3(co,[size(co,1), size(co,2)/20, size(co,3)/20]);
cross1=imresize3(cross, [size(cross,1), size(cross,2)/20, size(cross,3)/20]);
ref=single(sqrt(co1.^2+cross1.^2));
vol=ref(:,:,:);
sur=surprofile_agar(ref,'PSOCT',1); % find surface, downsample=20
Z_step=3;
sur_px = 5;
mus_depth = 150;

%Reshape data into 2D vector
start_depth = sur+sur_px; %*** fix start depth to be a surface 2D
fit_depth=min(size(vol,1)-5,start_depth+mus_depth-1);%***
vol2 = vol(start_depth:fit_depth,:,:);
aline_len = fit_depth-start_depth;
for i=1:size(vol2,2)
    for j=1:size(vol2,3)
        z_1D = transpose(((0:aline_len(i,j))+sur_px)*Z_step); % *** make sure that zdata are custom for each aline
        z_3D(:,i,j) = z_1D;%repmat(z_1D,1,N_alines);
    end
end

ydata = double(vol2(:));  % Replace with your actual ydata matrix
zdata = z_3D(:);  % Replace with your actual zdata matrix

% Initialize result
result = zeros(size(zdata));
aline_len = aline_len(:);
% Compute the result for each column in param4_matrix
zf_1D = zf_l(:);
for col = 1:length(zf_1D)
        pxS = 1 + (aline_len(col)+1)*(col-1);
        pxE = pxS + (aline_len(col));
        out = double(shift+sqrt(mub.* exp(-2 .* mus .* zdata(pxS:pxE)) .* ...
            (1 ./ (1 + ((zdata(pxS:pxE) - (zf_1D(col))) ./ zr).^2))));
        result(pxS:pxE) = out;
        y = ydata(pxS:pxE);
        r2(col) = 1 - sum((y - out).^2)/sum((y - mean(y)).^2);
end

figure(3); plot(result);hold on; plot(ydata);hold off; legend('fitting', 'data');

