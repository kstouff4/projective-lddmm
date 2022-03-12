function[] = getMRICoords_by_katie(subject,block)
%clear all;
close all;
fclose all;


addpath /cis/home/dtward/Functions/avwQuiet
addpath /cis/home/dtward/Functions/plotting
addpath /cis/home/dtward/Functions/dirn
addpath /cis/home/dtward/Functions/frame2Gif

% input names
%subject = 'Brain2';
%block = '1';
%block = '2';
% block = '3'; % most of my tuning was block 3
% NOTE: "brain3" is the same as "Brain2", preference is to use uppercase B version.
switch subject
    case 'Brain1'
        switch block
            case '1'
                mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/mri/Hip_Amyg_b0.img'];
                stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/histology/down_032/Brain1 AmyHippo/']; % Brain1 Hippo is alternative
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45];
            case '2'
                mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/mri/Hippocampus_b0_adjusted.img'];
                stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/histology/down_032/Brain1 Hippo/']; % Brain1 Hippo is alternative           
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43];

            otherwise
                disp('Defaulting to block 2')
                mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/mri/Hippocampus_b0_adjusted.img'];
                stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/histology/down_032/Brain1 Hippo/']; % Brain1 Hippo is alternative
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43];

        end
        stains = {'LFB H&E'};
        
        rigid_dir = '';
        sxy = 480;
        sz = 0.8*5;
        
    case 'Brain2'
        mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/mri/AD_Hip' block '/AD_Hip' block '_b0.img']; % should I do a bias correction first?
        stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/histology/down_032/AD_Hip' block '/'];
        stains = {'Amyloid','LFB','Tau'}; % these should be appended to the root dir
        sxy = 480;
        sz = 0.8;
        rigid_dir = '/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/brain3/histology/alignment_second_try/subject_brain3_block_1_initial_rigid_test0/';
        zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    case 'Brain3'
        switch block
            case '1'
                mri_filename = ['/cis/project/exvivohuman_11T/data/miller_mori_troncoso/20161206AD/20161206Amy/20161206Amy_b0.img'];
            case '2'
                mri_filename = ['/cis/project/exvivohuman_11T/data/miller_mori_troncoso/20161206AD/20161206HipAmy/20161206HipAmy_b0.img'];
            case '3'    
                mri_filename = ['/cis/project/exvivohuman_11T/data/miller_mori_troncoso/20161206AD/20161206Hip1/20161206Hip1_b0.img'];
            case '4'    
                mri_filename = ['/cis/project/exvivohuman_11T/data/miller_mori_troncoso/20161206AD/20161206Hip2/20161206Hip2-b0.img'];
        end
        % mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain3/aligned_block_' block '.img'];
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/histology/down_032/AD_Hip' block '/'];
        stains = {'Amyloid','Tau'};
        sxy = 480;
        sz = 0.8;
        rigid_dir = ''
        zInput = [1,2,3,4,5,6,7,8,9,10,11,12,13]; % should denote the spacing of the slices (i.e. intervals of 1 mm assumed between integers)
    case 'Brain5'
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/histology/down_032/AD_Hip' block '/'];
        stains = {'Amyloid','Tau'};
        sxy = 480;
        sz = 0.8;
        rigid_dir = '';
        tifstring = '*_crop_filterNS.tif';
        switch block
            case '1'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114amy-hippo_b0.img'];
            case '2'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114hippo1_b0.img'];
            case '3'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9]; % only use every 5 slices for which have amyloid and tau
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114hippo2_b0.img'];
                sz = 0.8*5;
                sz = 0.8;
        end
    
    otherwise
        disp('Defaulting to Brain 2')
        mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/mri/AD_Hip' block '/AD_Hip' block '_b0.img']; % should I do a bias correction first?
        stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/histology/down_032/AD_Hip' block '/'];
        stains = {'Amyloid','LFB','Tau'}; % these should be appended to the root dir
        sxy = 480;
        sz = 0.8;
        rigid_dir = '/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/brain3/histology/alignment_second_try/subject_brain3_block_1_initial_rigid_test0/';
        zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
end
down = 32;

output_prefix = ['subject_' subject '_block_' block '_alignment_diffeo_test9_by_katie_single/']; 

% start by loading MRI (assume 1 block only)
[avw,nxM,dxM,xM,yM,zM] = avw_img_read_domain([mri_filename]);
IM = avw.img;
% we need to zero pad; zero pad for all images
IM = padarray(IM,[1 1 1], 0, 'both');
nxM = nxM+2;
xM = [xM(1)-dxM(1), xM, xM(end)+dxM(1)];
yM = [yM(1)-dxM(2), yM, yM(end)+dxM(2)];
zM = [zM(1)-dxM(3), zM, zM(end)+dxM(3)];


IM = IM - mean(IM(:));
IM = IM / std(IM(:));

xM = xM - mean(xM);
yM = yM - mean(yM);
zM = zM - mean(zM);
[XM,YM,ZM] = meshgrid(xM,yM,zM);

indIm = IM>0;
xCent = (indIm).*XM;
xCent = sum(xCent(:))/sum(indIm(:));
yCent = (indIm).*YM;
yCent = sum(yCent(:))/sum(indIm(:));
zCent = (indIm).*ZM;
zCent = sum(zCent(:))/sum(indIm(:));
disp(num2str(xCent))
disp(num2str(yCent))
disp(num2str(zCent))

save([output_prefix 'coords.mat'],'xM','yM','zM','xCent','yCent','zCent');
end