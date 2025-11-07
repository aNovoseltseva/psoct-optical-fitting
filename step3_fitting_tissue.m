
%id='1'; %comment it when run fitting in parallel
% fitting scattering and birefringence of brain tissue
% run this after OCT_recon.m is finished
% clear all;close all;clc;
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

elseif res==3  % if tile size 3.3x3.3 mm use following parameters
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

% load parameters extracted from agar tile fitting (agar has constant mus and mub)
cd(strcat(folder,'fitting/'));
load('zf_zr_sim_fit.mat');

%% fit the entire volume

for islice=id%!!!!!!!
    cd(datapath);
    filename0=dir(strcat('co-',num2str(islice),'-*.dat')); 
    dirname=strcat(folder,'fitting/vol',num2str(islice));
    mkdir(dirname);
    for iFile=istart:istop
        name=strsplit(filename0(1).name,'.');  
        name_dat=strsplit(name{1},'-');
        slice_index=islice;
        % Xrpt and Yrpt are x and y scan repetition, default = 1
        Zsize = str2num(name_dat{4}); Xrpt = 1; Xsize=str2num(name_dat{5}); Yrpt = 1; Ysize = str2num(name_dat{6});
        dim1=[Zsize Xrpt Xsize Yrpt Ysize];     % tile size for reflectivity 
        name1=strcat(datapath,'co-',num2str(islice),'-',num2str(iFile),'-',num2str(Zsize),'-',num2str(Xsize),'-',num2str(Ysize),'.dat'); % gen file name for reflectivity
        if isfile(name1)
            % load reflectivity data
            co = ReadDat_int16(name1, dim1)./65535*4; 
            name1=strcat(datapath,'cross-',num2str(islice),'-',num2str(iFile),'-',num2str(Zsize),'-',num2str(Xsize),'-',num2str(Ysize),'.dat'); % gen file name for reflectivity
            cross = ReadDat_int16(name1, dim1)./65535*4; 
        else
            co=zeros(Zsize,Xsize,Ysize);
            cross=zeros(Zsize,Xsize,Ysize);
        end
        % pause here, use view3D(co) to get value for zf and
        message=strcat('Tile No. ',string(iFile),' is read.', datestr(now,'DD:HH:MM'),'\n');
        fprintf(message);
        Optical_fitting_AN_v4(co, cross, islice, iFile, folder,0, 130, 80, ds_factor, zf_l, zr, shift, 2);


    end

    metric_stitch('mus', P2path, folder,disp,mosaic,round(pxlsize./ds_factor),islice,pattern,sys,ds_factor,stitch);     
    metric_stitch('mub', P2path, folder,disp,mosaic,round(pxlsize./ds_factor),islice,pattern,sys,ds_factor,stitch);
    metric_stitch('bfg', P2path, folder,disp,mosaic,round(pxlsize./ds_factor),islice,pattern,sys,ds_factor,stitch);

    message=strcat('slice No. ',string(islice),' is fitted and stitched.', datestr(now,'DD:HH:MM'),'\n');
    fprintf(message);
end
system(strcat('chmod -R 777',{' '},folder));
clearvars -except id
