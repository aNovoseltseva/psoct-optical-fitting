%% ----------------------------------------- %%
% Note Mar 14 2023:
%
% Current version of code does FOV correction, grid correction &
% MIP/AIP/Retardance/Orientation generation, surface finding & profiling, 
% volume stitching
%
% Fitting is in Fitting_after_recon.m
%

%
% Note Jan 28 2020:
%
% Current version of code does FOV correction, grid correction & MIP/AIP
% generation, surface finding & profiling, mosacing & blending
%
% Suggested parallelization method: one slice per thread
% 
% volume stitching is currently in another code
%
% Note Oct 28 2019:
%
% Current version of code does FOV correction, grid correction & MIP/AIP
% generation, surface finding & profiling, mosacing & blending
% 
% volume stitching is currently in another code
%
% Note Nov 5 2019:
% 
% All parameters were moved to the beginning of the script
% volume stitching is currently in another code
%
% Note 09/02/2020
%
% Current version of code does FOV correction(optional), grid correction(optional) & MIP/AIP
% generation, surface finding & profiling, mosacing & blending
%
%
% Author: Jiarui Yang, Shuaibin Chang
%%%%%%%%%%%%%%%%%%%%%%%

%id='3';
%% parameters starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add subfunctions for the script. Change directory if not running on BU SCC
info = readlines('step0_info_d1.txt');

addpath(strcat(pwd, '\functions'));
% specify dataset directory
OCTpath  = info(2);  % OCT data path.              
P2path = info(5);  % 2P data path.
% leave unchanged if no 2P data is avalable                                                          

nslice=str2num(info(8));   % total number of imaging slices.                                                       
stitch=str2num(info(11)); % 1 means using OCT data to generate stitching coordinates, 
% 0 means using 2P stitching coordinates.                                                            ADJUST FOR EACH SAMPLE!!!
remove_agar = 0; % 1 means remove agar boundary in volume recon, 0 means keep agar boundary. 
% Use 1 when there is no depth scanning during imaging, i.e., image twice for 450um slice            ADJUST FOR EACH SAMPLE!!!
ten_ds = 0; % 1 means use 10x10 pixel for volume recon, which is 30x30x35 um voxel. 
% 0 means use 4x4 pixel, which is 12x12x14um voxel. For small sample size, use 0. For large sample size, use 1.
highres=1;% 1 means 2PM is imaged at 1um step size, 0 means 2PM is imaged at 2um step size          ADJUST FOR EACH SAMPLE!!!

% if using OCT images to generate stitching coordinates, you need three slices that are separated in z for stitching (why three? three channels in RGB format.)
% However, say, you have total of 100 slices, but you don't have enough space in SCC so you want to run 50 slices at a time, 
% then the following slice numbers should be smaller than 50
stitch_slice1=3;    % if stitch == 1, first slice for stitching                                      ADJUST FOR EACH SAMPLE!!!
stitch_slice2=12;   % if stitch == 1, second slice for stitching                                     ADJUST FOR EACH SAMPLE!!!
stitch_slice3=20;   % if stitch == 1, third slice for stitching                                      ADJUST FOR EACH SAMPLE!!!

njobs=1; % number of parallel tasks per slice in SCC parallel processing, default to be 1, meaning each task handles one slice
sys = 'PSOCT'; % specify OCT system name. default to be 'PSOCT
% specify mosaic parameters, you can get it from Imagej stitching
% xx yy is positive for dataset acquired after sep 06

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
xy=0;     % xy is the Y displacement of two adjacent tile align in the X direction, default to 0
yx=0;      % xx is the X displacement of two adjacent tile align in the Y direction, default to 0
numX=str2num(info(17));    % #tiles in X direction                                                                  
numY=str2num(info(20));   % #tiles in Y direction                                                                   

disp=[xx xy yy yx];
mosaic=[numX numY Xoverlap Yoverlap];

pattern = 'bidirectional';  % mosaic pattern, could be bidirectional or unidirectional

ntile=numX*numY;
%%%%%--------------------------ADJUST FOR EACH AMPLE!!!------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
start_pixel=48+15; % start depth for calculating MIP, retardance, and volume reconstruction. 
% Should be 10-20 pixels bellow tissue surface.                                                      

slice_thickness = 44; % imaging slice thickness in pixels. 
% For U01 sample should be 70, for Ann Mckee samples should be 44.                                   
aip_threshold=0.007; % intensity threshold for AIP (before BaSiC shading correction) 
% to remove agarous.                                                                                 
aip_threshold_post_BaSiC=aip_threshold/4; % intensity threshold for AIP (after BaSiC shading correction) to remove agarous. 

% the $SGE-TASK-ID environment variable read in is CHARACTER, need to transfer to number

id = str2num(id);                                                                                
section=ceil(ntile/njobs); % total tiles per each paralle task, usually equal to total tiles per slice
istart=1;% starting tile number for each parallel task         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
istop=section; % end tile number for each parallel task          !!!!!!!!!!!!!!  IF WANT TO REDO STITCHING WITHOUT REPEATING DISTORTION CORRECTION, SET ISTOP=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters stop here


% load distortion correction matrixes, Optional
if(strcmp(sys,'PSOCT')&&highres~=1)
    folder_distort='/projectnb/npbssmic/ns/distortion_correction/OCT_grid_after_20240418/';                     
    fileID = fopen(strcat(folder_distort,'grid matrix.bin'), 'r'); 
    grid_matrix = fread(fileID,'double');
    fclose(fileID);
    grid_matrix=reshape(grid_matrix, 4,pxlsize(1)+100,pxlsize(2)+100);
end
% directory that stores distortion corrected 3D tiles. Optional
corrected_path=strcat(OCTpath,'dist_corrected/'); 
mkdir(strcat(OCTpath,'dist_corrected'));
mkdir(strcat(OCTpath,'dist_corrected/volume'));
cd(OCTpath);

filename0=dir(strcat('*AB.dat')); % get file names   
% Read volume dimension before loading file
try
    name=strsplit(filename0(1).name,'.');  
    name_dat=strsplit(name{1},'-');
    % Xrpt and Yrpt are x and y scan repetition, default = 1
    Zsize = str2num(name_dat{3}); Xrpt = 1; Xsize=str2num(name_dat{4}); Yrpt = 1; Ysize = str2num(name_dat{5})/2;
    dim1=[Zsize Xrpt Xsize Yrpt Ysize*3];     % tile size for reflectivity 
    if strcmp(sys,'Thorlabs')
        dim1=[400 1 400 1 400];
    end
catch
    display('no raw data found')
    Zsize=300;
end
    
% create folder for AIPs and MIPs
create_dir(nslice, OCTpath);  
% function that finds tissue surface 


if ~isfile(strcat(OCTpath,'surface.mat'))
    %if you have only one tile 
    if ntile < 20
        % find curvature with 
        Find_FOV_Curvature(OCTpath);
    elseif id==3
        Hui_sum_all(OCTpath,ntile, 1,Zsize,Ysize); % parameters: data path, ntile, slice#, Z pixels. slice# should be the one that was cut even.
    else
        pause(7200)
    end
end
cd(OCTpath);
load(strcat(OCTpath,'surface.mat'));
if exist('curvature','var')
    surface = curvature;
end
surface=round(surface-min(surface(:)));
%%

for islice=id
     cd(OCTpath)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR LOOP TO PERFORM DISTORTION CORRECTION ON EACH TILE
    for iFile=istart:istop
        name1=strcat(num2str(islice),'-',num2str(iFile),'-',num2str(Zsize),'-',num2str(Xsize),'-',num2str(Ysize*2),'-AB.dat'); 
        
        folderPath = "/projectnb/npbssmic/ns/PMI_study_monkey_brain/PMI_Kiki/data/dist_corrected/";
        fileName = strcat('ori-',num2str(islice),'-',num2str(iFile),'-280-700-700.dat');
        fullPath = fullfile(folderPath, fileName);
    

        % load data
        ifilePath = strcat(OCTpath,name1);
        amp = ReadDat_int16(ifilePath, dim1)./65535*4;
        message=strcat('Tile No. ',string(iFile),' is read.', datestr(now,'DD:HH:MM'),'\n');
        fprintf(message);
        
        % distortion correction
        if(strcmp(sys,'PSOCT')) 
            % separating data into different channels
             cross=amp(:,101:end,1:Ysize); % neglect the first 100 pixel in x because of blade distortion
             co=amp(:,101:end,Ysize+1:Ysize*2);
             ori=amp(:,101:end,2*Ysize+1:Ysize*3);
             ori=ori.*90;
             ori(ori>180)=ori(ori>180)-180;
             % correct for tissue surface tilt
             co = FOV_curvature_correction(co, surface, size(co,1), size(co,2), size(co,3));
             aip=squeeze(mean(co(1:110,:,:),1)); %!!!!! stop here to check profile with view3D(co) !!!!!!!!!!!!!!!!!!!!!!!!
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Use the following line to find aip_threshold.
             % aip_threshold should be 0.01-0.02 larger than the average
             % pixel intensity of agarous
             
             % When running OCT_recon on SCC, DO COMMONT this line!!!
             % figure;imagesc(aip);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if (mean2(aip)>aip_threshold || std2(aip)>aip_threshold/4)
                 cross = FOV_curvature_correction(cross, surface, size(cross,1), size(cross,2), size(cross,3));  
                % correct distortion in X and Y dimension
                if highres~=1
                    cross = Grid_correction(cross, grid_matrix, 50+pxlsize(1), 51, 50+pxlsize(2), 51, size(cross,1));  
                    co = Grid_correction(co, grid_matrix, 50+pxlsize(1), 51, 50+pxlsize(2), 51, size(co,1));  
                end
                 %correct distortion in Z dimension 
                 ori = FOV_curvature_correction(ori, surface, size(ori,1), size(ori,2), size(ori,3)); 
                 if highres~=1
                    ori=ori(:,51:50+pxlsize(1),51:50+pxlsize(1));
                 end
             else
                 if highres~=1
                    cross = zeros(size(co,1),pxlsize(1),pxlsize(2));  
                    co=co(:,51:50+pxlsize(1),51:50+pxlsize(2)); 
                 end
                 ori=zeros(size(co,1),pxlsize(1),pxlsize(2)); 
             end
             ref=sqrt(cross.^2+co.^2);
             ret=atan(cross./co)./pi*180;
        end
         % surface profiling and save to folder 
         sur=surprofile2(ref,sys,10); % to check how well blade cuts                          
         surname=strcat(OCTpath,'surf/vol',num2str(islice),'/',num2str(iFile),'.mat');
         save(surname,'sur');
         % Generating AIP, MIP, retardance AIP in mat
         %mip -max intentsity proejction
         %aip- average intensity projection
         mip=squeeze(max(ref(start_pixel:start_pixel+50,:,:),[],1)); 
         aip=squeeze(mean(ref(1:110,:,:),1));
         ret_aip=squeeze(mean(ret(start_pixel:start_pixel+100,:,:),1));  
         
         [X Y]=meshgrid(0.1:pxlsize(1)/10); 
         if highres==1
            [Xv Yv]=meshgrid(0.1:0.0992:69.5);
         else
            [Xv Yv]=meshgrid(0.1:0.0991:99.1);
         end
         Vq=interp2(X,Y,sur,Xv,Yv);
         sur=round(Vq);
         ori2D=Gen_ori_2D(ori,sur,40);
         
         
         % saving distortion-corrected tiles. Optional not saving agr
         %if (mean2(aip)>aip_threshold || std2(aip)>aip_threshold/4)
           name1=strcat(num2str(islice),'-',num2str(iFile),'-',num2str(size(cross,1)),'-',num2str(size(cross,2)),'-',num2str(size(cross,3)),'.dat'); % gen file name for reflectivity
           FILE_ref=strcat(corrected_path, 'cross-', name1);
           FID=fopen(FILE_ref,'w');
           fwrite(FID,cross./4*65535,'uint16');
           fclose(FID);

           name2=strcat(num2str(islice),'-',num2str(iFile),'-',num2str(size(co,1)),'-',num2str(size(co,2)),'-',num2str(size(co,3)),'.dat'); % gen file name for reflectivity
           FILE_ret=strcat(corrected_path, 'co-', name2);
           FID=fopen(FILE_ret,'w');
           fwrite(FID,co./4*65535,'uint16');
           fclose(FID);

           name3=strcat(num2str(islice),'-',num2str(iFile),'-',num2str(size(ori,1)),'-',num2str(size(ori,2)),'-',num2str(size(ori,3)),'.dat'); % gen file name for reflectivity
           FILE_ret=strcat(corrected_path, 'ori-', name2);
           FID=fopen(FILE_ret,'w');
           fwrite(FID,(ori./180).*65535,'uint16');
           fclose(FID);
         %end
         
        % Saving AIP.tif
        tiffname=strcat(OCTpath,'aip/vol',num2str(islice),'/','AIP.tif');
        SaveTiff(aip,iFile,tiffname);
        
       % Saving retardance AIP.tif
        tiffname=strcat(OCTpath,'retardance/vol',num2str(islice),'/','RET.tif');
        SaveTiff(ret_aip,iFile,tiffname);

        % Saving orientation ORI.tif
        tiffname=strcat(OCTpath,'orientation/vol',num2str(islice),'/','ORI.tif');
        SaveTiff(ori2D,iFile,tiffname);

        % Saving MIP.tif
        tiffname=strcat(OCTpath,'mip/vol',num2str(islice),'/','MIP.tif');
        SaveTiff(mip,iFile,tiffname);

        fprintf(strcat('Tile No. ',string(iFile),' is reconstructed.', datestr(now,'DD:HH:MM'),'\n'));

    end  
    % log file confirming this slice has finished distortion correction
    fid=fopen(strcat(OCTpath,'dd',num2str(id),'.txt'),'w');
    fclose(fid);
    system(strcat('chmod -R 777', {' '},OCTpath));

%% Stitching
    % pause untile the stitching slices finish distortion correction
    if id==1 && stitch==1 && ~isfile(strcat(OCTpath,'aip/RGB/TileConfiguration.registered.txt'))
        display('wait for distortion recon to be finished')
        while ~(isfile(strcat(OCTpath,'dd',num2str(stitch_slice1),'.txt')) && isfile(strcat(OCTpath,'dd',num2str(stitch_slice2),'.txt')) && isfile(strcat(OCTpath,'dd',num2str(stitch_slice3),'.txt')))
            pause(600);
        end
        %if we dont have 2P data use this code for stitching
        display('generating stitching coordinates')
        Gen_OCT_coord(OCTpath,disp,mosaic,pattern, stitch_slice1, stitch_slice2, stitch_slice3); % generating OCT stitching coordinates
    end
    % pause until stitching cooridates are generated
    if stitch==1
        display('wait for stitching coordinates to be generated')
        while ~isfile(strcat(OCTpath,'aip/RGB/TileConfiguration.registered.txt')) 
            pause(600);
        end
    end
    display('start stitching')
    %AIP_stitch(P2path,OCTpath,disp,mosaic,pxlsize,islice,pattern,sys,stitch,aip_threshold,aip_threshold_post_BaSiC, highres);               % stitch AIP
    %MIP_stitch('mip', P2path,OCTpath,disp,mosaic,pxlsize,islice,pattern,sys,stitch,aip_threshold_post_BaSiC, highres);                      % stitch MIP
    %Surf_stitch('sur',P2path,OCTpath,disp,mosaic,pxlsize/10,islice,pattern,sys,stitch, highres);                                            % stitch surface
    %Ret_stitch('ret_aip', P2path,OCTpath,disp,mosaic,pxlsize,islice,pattern,sys,stitch,aip_threshold_post_BaSiC, highres);                  % stitch retardance AIP
    
    system(strcat('chmod -R 777', {' '},OCTpath));
    %Ori_stitch('ori2D', P2path,OCTpath,disp,mosaic,pxlsize,islice,pattern,sys,stitch,aip_threshold_post_BaSiC, highres);                    % stitch orientation
    %Gen_ori_RGB(OCTpath,islice, 0.015);
    % this is for stitching volume

    BaSiC_shading_and_ref_stitch(islice,P2path,OCTpath, numX*numY, start_pixel, slice_thickness,stitch, Xoverlap, Yoverlap, highres);  % volume recon 
    fprintf(strcat('Slice No. ',num2str(islice),' is stitched.', datestr(now,'DD:HH:MM'),'\n'));



end
%% log file indicates this slice has done reconstruction
fid=fopen(strcat(OCTpath,'dist_corrected/volume/log',num2str(id),'.txt'),'w');
fclose(fid);
cd(strcat(OCTpath,'dist_corrected/volume/'))
logfiles=dir(strcat(OCTpath,'dist_corrected/volume/log*.txt')); 
% The last slice finish volume recon will do stacking 
if length(logfiles)==nslice
    delete log*.txt
    Concat_ref_vol(nslice,OCTpath, aip_threshold_post_BaSiC, remove_agar, ten_ds);
    if ten_ds == 0
        depth_normalize_2(OCTpath, 50,50, 0); % volume intensity correction
    else
        depth_normalize_2(OCTpath, 20,20, 1); % volume intensity correction
    end
end
system(strcat('chmod -R 777',{' '},OCTpath));