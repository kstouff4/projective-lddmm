%{

Katie Notes:
1) This program aims to align scatter transformed histology images deformed
to 2D slices via rigid + nonrigid 3D transformations

 Example for running: 

%}
% 
%
% run align_v01_diffeo_by_katie_single('Brain1','1','Tau',1000,'_fAI_F_020222.tif','020222')
%

function[] = align_v01_diffeo_by_katie_single_round2(subject,block,stainC,iterat,suff,dateS)
clearvars -except subject block iterat lab stainC suff dateS;
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
% Note: zInput maps the first slice to 1 but skips the corresponding other
% slices
switch subject
    case 'Brain1'
        switch block
            case '1'
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain1/MRI/Hip_Amyg_b0.img'];
                stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain1/histology/down_032/AD_Hip1/']; % Brain1 Hippo is alternative
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45];
            case '2'
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain1/MRI/Hippocampus_b0_adjusted.img'];
                stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain1/histology/down_032/AD_Hip2/']; % Brain1 Hippo is alternative           
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43];

            otherwise
                disp('Defaulting to block 2')
                mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/mri/Hippocampus_b0_adjusted.img'];
                stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain1/histology/down_032/Brain1 Hippo/']; % Brain1 Hippo is alternative
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43];

        end
        stains = {stainC};
        
        rigid_dir = '';
        sxy = 480;
        sz = 0.8*5; % approximately 0.2 microns apart (this amounts to 0.25)
        tifstring = '*.tif';
        
    case 'Brain2'
        mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/MRI/AD_Hip' block '_b0.img']; % should I do a bias correction first?
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain2/histology/down_032/AD_Hip' block '/'];
        stains = {stainC}; % these should be appended to the root dir
        sxy = 480;
        sz = 0.8;
        rigid_dir = '/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/brain3/histology/alignment_second_try/subject_brain3_block_1_initial_rigid_test0/';
        switch block
            case '1'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
            case '2'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
            case '3'
            zInput = [1, 2, 3, 4, 5, 6, 7];
        end
        tifstring = ['*' suff];
    
    case 'Brain3'
        switch block
            case '1'
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/MRI/20161206Amy_b0.img'];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
            case '2'
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/MRI/20161206HipAmy_b0.img'];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
            case '3'    
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/MRI/20161206Hip1_b0.img'];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
            case '4'    
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/MRI/20161206Hip2-b0.img'];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
        end
        %mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain3/aligned_block_' block '.img'];
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain3/histology/down_032/AD_Hip' block '/'];
        stains = {stainC};
        sxy = 480;
        sz = 0.8;
        rigid_dir = '';
        tifstring = '*_crop_filterNS.tif';
        
    case 'Brain4'
        mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/MRI/aligned_block_' block '.img'];
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/histology/down_032/AD_Hip' block '/'];
        switch block
            case '1'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/MRI/20161122Hy-Amy_b0.img'];
            case '2'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/MRI/20161122Hip1_b0.img'];
            case '3'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain4/MRI/20161122Hip2_b0.img'];
        end
        stains = {stainC};
        sxy = 480;
        sz = 0.8;
        rigid_dir = '';
        tifstring = ['*' suff];
    case 'Brain5'
        stain_root_dir = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/histology/down_032/AD_Hip' block '/'];
        stains = {stainC};
        sxy = 480;
        sz = 0.8;
        rigid_dir = '';
        tifstring = ['*' suff];
        switch block
            case '1'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114amy-hippo_b0.img'];
            case '2'
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114hippo1_b0.img'];
            case '3'
                stains = {stainC}; % when trying other slices, don't use amyloid for alignment necessarily
                zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40];
                zInput = [1, 2, 3, 4, 5, 6, 7, 8]; % only use every 5 slices for which have amyloid and tau (3,7,12,17,22,27,32,37)
                mri_filename = ['/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/MRI/20161114hippo2_b0.img'];
                sz = 0.8*5;
                sz = 0.8;
        end

    otherwise
        disp('Defaulting to Brain 2')
        mri_filename = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/mri/AD_Hip' block '_b0.img']; % should I do a bias correction first?
        stain_root_dir = ['/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/Brain2/histology/down_032/AD_Hip' block '/'];
        stains = {stainC}; % these should be appended to the root dir
        sxy = 480;
        sz = 0.8;
        rigid_dir = '/cis/home/dtward/Documents/exvivohuman_11T/more_blocks/brain3/histology/alignment_second_try/subject_brain3_block_1_initial_rigid_test0/';
        zInput = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        tifstring = ['*' suff];
end
down = 32;


% for output
% repeat this for katie, save outputs
output_prefix = ['subject_' subject '_block_' block '_alignment_diffeo_test9_by_katie_single_round2/' num2str(iterat) '_iters_' dateS '/'];  
% katie's slices were generated with an output prefix on this line.

[mydir,myname,myext] = fileparts(output_prefix);
if ~exist(mydir,'dir')
    mkdir(mydir)
end

frame = [];
atlasframe = [];
framefignum = 99;
paramfignum = 98;
atlasfignum = 97;




%%
% my parameters
% optimization iterations
niter = 400;
% matching cost
sigmaM = 1.0;
% for diffeo
nt = 3;

% optimization params
% for 2d rigid
epTJ = 5e-1; % note 1e0 seems good but a bit on the big side
epRJ = 2e-8;
% I'm pretty happy with the above numbers
% trying 0 to test some things
% epTJ = 0;
% epRJ = 0;

% xy scale
epSxy = 5e-10;
epSxy = 2e-10;
epSxy = 1e-10; % reduction by 2 by Katie

% z scale
epSz = 1e-8;

% this z scale change is quite slow
epSz = 5e-8; % this might be too big
epSz = 2e-8; % still to big, start cutting off after 10
epSz = 1e-8;
% test 3
epSz = 5e-9;
% test 4
epSz = 2e-8;
% test 5
% epSz = 0;
% later tests 6 7
epSz = 2e-8;
% I want to make it smaller, but do different restarts
epSz = 1e-8;
epSz = 5e-9;
epSz = 2.5e-9; % reduction by 2 by Katie

% epSz = 0;
% epSxy = 0;

% 3d rigid
epR = 1e-9; % e-9 or e-10 are probably appropriate
epR = 0.5e-9; % reduction by 2 by Katie

% when I set epRJ to 0, I see that this rotation is actually quite big
epT = 1e-7; % was oscilating with 1e-7
epT = 1e-6; % still too big
epT = 1e-8; % oops looks like I accidentally made it bigger instead of smaller
epT = 0.5e-8; % reduction by 2 by Katie

% differential operators
p = 2;
a = 0.5; % half mm
% try twice as big
a = 1.0; % 1mm

% diffeo
sigmaR = 1e3;
% got things kind of working with basically no reg, now try to bring in reg
sigmaR = 1e2;
% still nothing coming up in reg
sigmaR = 1e1;
%

% update diffeo
epV = 1e-6;
epV = 1e-3; % try to make it too big
% still too low
epV = 1e-1;
% still too low
epV = 1e1;
% above is too much but just a bit too much
epV = 1e0;
% we kept getting oscilatoin
% try this, if it doesn't work, then make all parameters smaler
epV = 5e-1;
% it didn't work
epV = 1e0;
epV = 5e-1;  % i included a step_mult here of 0.5 that I had previously forgot
% note detA was about 2e5
epV = 2e-6;
epV = 1e-6; % reduction by 2 by Katie

if strcmp(block,'2')
    epSz = 2e-9;
    epSz = 1e-9;
    epSz = 0.5e-9; % reduction by 2 by Katie
    niter = 1000;
end
if strcmp(block,'3')
    epSz = 2e-9;
    epSz = 1e-9;
    epSz = 0.5e-9; % reduction by 2 by Katie
    niter = 1000;
%     epTJ = 2e-1; % this needed to be made smaller to stop oscilatoin
    
    % note I'm stlil getting oscilations now with diffeo, 
    % these are NOT in the reg energy, they are in the other parameters
    % maybe reducing all step sizes would be a good idea?
    
    % okay I still oscillated, but it was all due to LFB
    % I reduced lfb size and will leave everything else the same
    % by iteration 50 I was oscilating, it was the rotation angle
    % oscilating mostly
    % but again it was just one purple slice
    % so moving purple down again to 0.25;
    
    % still not working well
    % I think the problem may be the slice thickness again
    % let's give it another try 
    epSz = 5e-10;
    epSz = 2.5e-10; % reduction by 2 by Katie
    
    % better but still getting oscilations after a while
    
end
% let's just use those parameters for everyone
epSz = 2e-9;
epSz = 1e-9;
niter = 1000;
epSz = 5e-10;
epSz = 2.5e-10; % reduction by 2 by Katie
% these parameters start oscilating at around iter 200-300
% seems better now, one I reduced multiplier
% and two I padded the mri
niter = 3000;
% increase diffeo step
epV = 5e-6;
epV = 2.5e-6; % reduction by 2 by Katie
% increase reg
sigmaR = 1e0;
% increase reg
sigmaR = 1e-1;
% and increase epV
epV = 1e-5;
epV = 0.5e-5; %reduction by 2 by Katie
% decrease sigmaR
sigmaR = 1e-2;
% the above converged and didn't really go anywhere
sigmaR = 2e-2; % version 8
sigmaR = 5e-2; % version 9

% multiply lfb step size in slice matching
lfb_factor = 0.5;
lfb_factor = 0.333;
lfb_factor = 0.25;

% global stepsize multiplier
% for linear intensity transform
step_mult = 1.0;
% for cubic
step_mult = 0.5;
% % try reducing everything (haven't tried this yet)
% % step_mult = 0.25;
if (subject == 'Brain5')
    step_mult = 0.1; % try for no contrast transform
    sigmaR = 1e-2;
end
step_mult_min = step_mult/10;
% let's try a little less (test8)
step_mult_min = step_mult/5;

% more iterations
niter = 5000;
niter = iterat;

% rigid dir (empty string for do not load)
%rigid_dir = ''; % Brain1 parameter
%%
% start by loading MRI (assume 1 block only)
[avw,nxM,dxM,xM,yM,zM] = avw_img_read_domain([mri_filename]); % xM,yM, and zM already in mm depending on file
IM = avw.img;
% we need to zero pad; zero pad for all images
IM = padarray(IM,[1 1 1], 0, 'both');
nxM = nxM+2;
xM = [xM(1)-dxM(1), xM, xM(end)+dxM(1)];
yM = [yM(1)-dxM(2), yM, yM(end)+dxM(2)];
zM = [zM(1)-dxM(3), zM, zM(end)+dxM(3)];

disp('MRI values');
disp(mean(IM(:)));
disp(max(IM(:)));
disp(min(IM(:)));


% keep the same image values as before 
IMorig = IM - mean(IM(:));
IMorig = IMorig / std(IMorig(:));

IM = IMorig - min(IMorig(:));
IM = IM / (max(IMorig(:)) - min(IMorig(:)));
IM = IMorig;

xM = xM - mean(xM);
yM = yM - mean(yM);
zM = zM - mean(zM);
[XM,YM,ZM] = meshgrid(xM,yM,zM);
% danfigure(1);
% sliceView(xM,yM,zM,avw.img)

% compute the gradient
[IM_x,IM_y,IM_z] = gradient(IM,dxM(1),dxM(2),dxM(3));

vtx = zeros([size(IM),nt]); % maybe different init, but need size so can't do right at beginning
vty = zeros([size(IM),nt]);
vtz = zeros([size(IM),nt]);


%%
% now load histology
% danfigure(2);
dxJ = [down,down]; % difference in both directions is 32 pixels

themin = 1e10*[1,1];
themax = -1e10*[1,1];
max_X = -1e10;
max_Y = -1e10;
for s = 1 : length(stains)
    stain = stains{s};
    mydir = [stain_root_dir stain '/'];
    files = dirn([mydir tifstring]);
    n = length(files);
    disp('number of files')
    disp(n)
    for i = 1 : n
        tmp = imread([mydir files(i).name]); % rows by columns
%         imagesc(tmp)
        disp(size(tmp));
        %J{s,i} = (double(tmp)-min(tmp(:)))/(max(tmp(:)) - min(tmp(:))); % rescale to 0 to 1
        J{s,i} = double(tmp)-mean(tmp(:));
        J{s,i} = J{s,i}/std(J{s,i}(:)); % scale to pos and neg like MRI
        nxJ{s,i} = [size(tmp,2),size(tmp,1)]; % columns by rows 
        zJ(i) = zInput(i); % Assume same for all stains 
        xJ{s,i} = (0 : nxJ{s,i}(1)-1)*dxJ(1); % x coords (for intensity vals J)
        yJ{s,i} = (0 : nxJ{s,i}(2)-1)*dxJ(2); % y coords (for intensity vals J)
        max_X = max(max(nxJ{s,i}(1)),max_X); % indicates maximum size
        max_Y = max(max(nxJ{s,i}(2)),max_Y);
        xJ{s,i} = xJ{s,i} - mean(xJ{s,i});
        yJ{s,i} = yJ{s,i} - mean(yJ{s,i});
        [XJ{s,i},YJ{s,i}] = meshgrid(xJ{s,i},yJ{s,i});
        
        themin = min([themin;[min(xJ{s,i}),min(yJ{s,i})]]);
        themax = max([themin;[max(xJ{s,i}),max(yJ{s,i})]]);
        
        % initialize transform
        RJ{s,i} = eye(2);
        TJ{s,i} = [0;0];
        
        % initialize variables to save
        Bsave{i} = eye(4);
    end
end
zJ = zJ - mean(zJ);
% maximal domain that I will align all my slices within
xJ0 = themin(1) : dxJ(1) : themax(1);
yJ0 = themin(2) : dxJ(2) : themax(2);
zJ0 = zJ;
[XJ0,YJ0,ZJ0] = meshgrid(xJ0,yJ0,zJ0); % where to sample final MRI slices in?
xJ0_sing = (0 : max_X + 10 - 1)*dxJ(1);
yJ0_sing = (0 : max_Y + 10 - 1)*dxJ(2);
xJ0_sing = xJ0_sing - mean(xJ0_sing);
yJ0_sing = yJ0_sing - mean(yJ0_sing);
[XJ0_sing,YJ0_sing] = meshgrid(xJ0_sing,yJ0_sing); % x and y to use whereas zJ0 should be used 

% load rigid?
if ~isempty(rigid_dir)
    load([rigid_dir 'RJ.mat'])
    load([rigid_dir 'TJ.mat'])
end


%%
% differential operators
[FX,FY,FZ] = meshgrid((0:nxM(1)-1)/nxM(1)/dxM(1),(0:nxM(2)-1)/nxM(2)/dxM(2),(0:nxM(3)-1)/nxM(3)/dxM(3));
LL = ( 1.0 - a^2*2*( (cos(2.0*pi*FX*dxM(1))-1.0)/dxM(1)^2 + (cos(2.0*pi*FY*dxM(2))-1.0)/dxM(2)^2 + (cos(2.0*pi*FZ*dxM(3))-1.0)/dxM(3)^2  ) ).^(2*p);
K = 1.0./LL;


%%
% start looping
Esave = [];
EMsave = [];
sxysave = [];
szsave = [];
tsave = [];
thetasave = [];
ERsave = [];

% which frames to save? I want something log spaced
% don't want to save them all, just too big

nframes = 100; % KATIE CHANGE from 200
nlog = nframes;
while 1
    frames_to_save = unique(round(logspace(log(1)/log(10),log(niter)/log(10),nlog)));
    if length(frames_to_save) == nframes
        break
    else
        nlog = nlog + 1;
    end
end

for it = 1 : niter
    %%
    % the forward model
    % it will have several sampling stages
    % first diffeomorphism applied to atlas space
    % we will start with identity
    % then rigid transform will be sampled in slice space
    % then 2D rigid rescale sampled in tif space
    
    if mod(it,10) == 0
        disp(it)
    end
    % start with diffeo
    dt = 1/nt;
    phiinvx = XM;
    phiinvy = YM;
    phiinvz = ZM;
    IMt = repmat(IM,[1,1,1,nt]);
    for t = 1 : nt
        % deform image
        if t > 1
            F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
            IMt(:,:,:,t) = F(phiinvy,phiinvx,phiinvz);
        end
        % update phi
        Xs = XM - vtx(:,:,:,t)*dt;
        Ys = YM - vty(:,:,:,t)*dt;
        Zs = ZM - vtz(:,:,:,t)*dt;
        % subtract and add identity
        F = griddedInterpolant({yM,xM,zM},phiinvx-XM,'linear','nearest');
        phiinvx = F(Ys,Xs,Zs) + Xs;
        F = griddedInterpolant({yM,xM,zM},phiinvy-YM,'linear','nearest');
        phiinvy = F(Ys,Xs,Zs) + Ys;
        F = griddedInterpolant({yM,xM,zM},phiinvz-ZM,'linear','nearest');
        phiinvz = F(Ys,Xs,Zs) + Zs;
    end
    % cost
    vtxhat = fft(fft(fft(vtx,[],1),[],2),[],3);
    vtyhat = fft(fft(fft(vty,[],1),[],2),[],3);
    vtzhat = fft(fft(fft(vtz,[],1),[],2),[],3);
    ER = sum(sum(sum(sum( bsxfun(@times, (abs(vtxhat).^2 + abs(vtyhat).^2 + abs(vtzhat).^2) ,  LL ) ))))*prod(dxM)*dt/numel(vtx(:,:,:,1));
    ERsave = [ERsave,ER];
    
    F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
    phiIM = F(phiinvy,phiinvx,phiinvz);
    danfigure(atlasfignum);
    sliceView(xM,yM,zM,phiIM) % this will be important to visualize, also my data in this space
    if any(it==frames_to_save)
        atlasframe = [atlasframe, getframe(atlasfignum)];
        frame2Gif(atlasframe,[output_prefix 'deformed_atlas.gif'])
    end
    
    if mod(it,250) == 0
        phisaveX = phiinvx;
        phisaveY = phiinvy;
        phisaveZ = phiinvz;
        save([output_prefix 'relTrans_temp.mat'],'Bsave','phisaveX','phisaveY','phisaveZ');
    end
        
    
    %%
    % next stage is rigid to slice orientation
    % likely we will want to calculate gradients previously and transform the
    % gradients as well
    if it == 1
        RM = eye(3);
        RM = [1,0,0;
            0,0,1;
            0,-1,0];
        TM = [0;0;0];
        if strcmp(block, '1')
            switch subject
                case 'Brain1'
                    % Brain1
                    RM = [0, -1, 0;
                        1, 0, 0;
                        0, 0, 1];
                    TM = [0;0;0];
                case 'Brain2'
                    RM = [0, -1, 0;
                        0, 0, -1;
                        -1, 0, 0];
                    TM = [0;0;1];
                case 'Brain3'
                    RM = [0, -1, 0;
                        0,0,1;
                        -1,0,0];
                    TM = [0;0;0];
                case 'Brain5'
                    RM = [0, -1, 0;
                        -1,0,0;
                        0,0,-1];
                    TM = [0;0;0];
                otherwise
                    % Default to Brain2
                    RM = [0, -1, 0;
                        0, 0, -1;
                        -1, 0, 0];
                    TM = [0;0;1];
            end
            % Leave translations the same for each brain
%             TM = [0;0;0.5]; % bigger number makes you slide from left to right in my array
            
            % just see if a nice initialization will help avoid the z scale
            % exploding down at the beginning
%             TM = [-0.38998 -0.078941 1.0893]';
%             RM = [
%                 0.0975   -0.9859    0.1355
%                 -0.0708   -0.1426   -0.9873
%                 -0.9927   -0.0868    0.0837
%                 ];
            
        elseif strcmp(block,'2') 
            TM = [0;0;0];
            switch subject
                case 'Brain1'
                    % Brain 1
                    RM = [-1, 0, 0;
                        0, 1, 0;
                        0, 0, 1];
                case 'Brain2'
                    % Brain 2
                    RM = [0,0,-1;
                        -1,0,0;
                        0,1,0];
                case 'Brain3'
                    RM = [0,1,0;
                        0,0,1;
                        1,0,0];
                case 'Brain5'
                    RM = [0,1,0;
                        0,0,-1;
                        -1,0,0];
                otherwise
                    % Brain 2
                    RM = [0,0,-1;
                        -1,0,0;
                        0,1,0];
            end

        elseif strcmp(block,'3')
            TM = [0;0;0];
            switch subject
                case 'Brain2'
                     RM = [0,1,0;
                        1,0,0;
                        0,0,-1];
                     TM = [-1;-1;5];
                case 'Brain3'
                    RM = [0,-1,0;
                        0,0,1;
                        -1,0,0];
                case 'Brain5'
                    RM = [0,-1,0;
                        0,0,1;
                        -1,0,0];
                otherwise
                     RM = [0,1,0;
                        1,0,0;
                        0,0,-1];
            end
        elseif strcmp(block,'4')
            switch subject
                case 'Brain3'
                    RM = [0,1,0;
                        0,0,1;
                        1,0,0];
                otherwise
                    RM = [0,1,0;
                        0,0,1;
                        1,0,0];
            end
            TM = [0;0;0];


        end
    end
    AR = [RM,TM;0,0,0,1]; % R for rigid
    BR = inv(AR);
      
    ASz = diag([1,1,sz,1]);
    
    % now xy
    ASxy = diag([sxy,sxy,1,1]);
    
    
    
    
    %%
    % now map to exact images, and calculate cost
    TJbar = [0;0];
    thetaJbar = 0;
    for s = 1 : length(stains)
        for i = 1 : n
            
            % map from R3 to R2,0 so it can be invertible
            %AJ = [[RJ{s,i}, [0;0];0,0,1], [TJ{s,i};0];[0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]]; % no 2D transformations here
            BJ = inv(AJ);
            danfigure(6);
            
            subplot(2,2,1)
            imagesc(xJ{s,i},yJ{s,i},J{s,i})
            axis image
            title('Observed image')
            
            % I could keep pushing back
            A = AJ * ASxy * ASz * AR;
            B = inv(A);
            
            Xs = B(1,1)*XJ{s,i} + B(1,2)*YJ{s,i} + B(1,3)*zJ(i) + B(1,4); % a1*x + b1*y + c1*z
            Ys = B(2,1)*XJ{s,i} + B(2,2)*YJ{s,i} + B(2,3)*zJ(i) + B(2,4);
            Zs = B(3,1)*XJ{s,i} + B(3,2)*YJ{s,i} + B(3,3)*zJ(i) + B(3,4);
            
            F = griddedInterpolant({yM,xM,zM},phiinvx-XM,'linear','nearest');
            phiinvBx = F(Ys,Xs,Zs) + Xs;
            F = griddedInterpolant({yM,xM,zM},phiinvy-YM,'linear','nearest');
            phiinvBy = F(Ys,Xs,Zs) + Ys;
            F = griddedInterpolant({yM,xM,zM},phiinvz-ZM,'linear','nearest');
            phiinvBz = F(Ys,Xs,Zs) + Zs;
            
            F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
            IS{s,i} = F(phiinvBy,phiinvBx,phiinvBz);
            
            %AJs =[[eye(2), [0;0];0,0,1], [TJ{s,i};0];[0,0,0,1]];
            AJs = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];
            As = AJs * ASxy * ASz * AR;
            Bs = inv(As);
            Xss = Bs(1,1)*XJ{s,i} + Bs(1,2)*YJ{s,i} + Bs(1,3)*zJ(i) + Bs(1,4); % a1*x + b1*y + c1*z
            Yss = Bs(2,1)*XJ{s,i} + Bs(2,2)*YJ{s,i} + Bs(2,3)*zJ(i) + Bs(2,4);
            Zss = Bs(3,1)*XJ{s,i} + Bs(3,2)*YJ{s,i} + Bs(3,3)*zJ(i) + Bs(3,4);
            phisaveXb{s,i} = F(Yss,Xss,Zss) + Xss;
            phisaveYb{s,i} = F(Yss,Xss,Zss) + Yss;
            phisaveZb{s,i} = F(Yss,Xss,Zss) + Zss;
            % I also  may want to consider a slice specific diffeo at this
            % point
            
            subplot(2,2,2)
            imagesc(xJ{s,i},yJ{s,i},IS{s,i})
            
            axis image
            title('Aligned MRI')
            colormap gray
             
            EMsi(s,i) = sum(sum(sum((IS{s,i} - J{s,i}).^2)))*prod(dxJ)/2.0/sigmaM^2;
            %                 disp(EMsi(s,i))
            %                 return
            
            %             disp(EMsi(s,i))
            %     end
            % end
            % disp(EM)
            
            % %%
            % now the forward model is finished
            % the first optimization is already done (the unkown intensit ymodel)
            epTJ_ = epTJ;
            epRJ_ = epRJ;
            if s == 2 % this stain has higher gradients                
                epTJ_ = epTJ*lfb_factor;
                epRJ_ = epRJ*lfb_factor;
            end
            
            %             danfigure(7);
            % for s = 1 : length(stains)
            %     for i = 1 : n
            
            err_rgb = IS{s,i} - J{s,i};
            %             subplot(2,2,1)
            %             imagesc(err+0.5) % maximum of error is 0.5, minimum is -0.5
            %             axis image
            %
            
            
            
            
            
            
            % next stage
            %tmp = reshape( DDfun(IS{s,i}(:)) * coeffs{s,i} , size(fIS{s,i}));
            %err = sum(err_rgb .* tmp, 3); % this sums over rgb
            err = err_rgb;
            %             danfigure(7);
            %             subplot(2,2,2)
            %             imagesc(err) % maximum of error is 0.5, minimum is -0.5
            %             axis image
            
            % image gradient
            [I_x,I_y] = gradient(IS{s,i},dxJ(1),dxJ(2));
            
            % NOT UPDATING ON THIS SECOND ROUND
            % get gradients
            %AJ = [RJ{s,i},TJ{s,i};[0,0,1]];
            %rgrad = zeros(2,3);
            %{
            for i0 = 1 : 2
                for i1 = 1 : 3
                    %                 deltaR = [0,0,1;0,0,0;0,0,0];
                    deltaR = [double((1:2)==i0)' * double((1:3)==i1);[0,0,0]];
                    RDeltaRRinv = AJ * (deltaR /AJ);
                    Xs = RDeltaRRinv(1,1)*XJ{s,i} + RDeltaRRinv(1,2)*YJ{s,i} + RDeltaRRinv(1,3);
                    Ys = RDeltaRRinv(2,1)*XJ{s,i} + RDeltaRRinv(2,2)*YJ{s,i} + RDeltaRRinv(2,3);
                    rgrad(i0,i1) = -sum(sum(  err .* (I_x.*Xs + I_y.*Ys)  ))*prod(dxJ);
                end
            end
            rotgrad = rgrad(1:2,1:2) - rgrad(1:2,1:2)';
            tgrad = rgrad(1:2,3);
            
            % looks like either method works
            %         TJ{s,i} = TJ{s,i} - epTJ*tgrad;
            %         RJ{s,i} = RJ{s,i} * expm(-epRJ*rotgrad); % right perturbation
            AJ_right_grad_mult{s,i} = expm([[-step_mult*epRJ_*rotgrad,-step_mult*epTJ_*tgrad];[0,0,0]]);
            
            
            danfigure(framefignum);
            %                 set(framefignum,'position',[100 100 1024*2 512*2]);
            set(framefignum,'position',[100 100 1024*2 768]);
%             set(framefignum,'position',[100 100 1024 768/2]); % too many frames so I made it smaller, now I'm using less frames
            subplotdan(length(stains)*3,n,(s-1)*n*3 + i);
            imagesc(J{s,i})
            axis image
            axis off
            subplotdan(length(stains)*3,n,(s-1)*n*3 + n + i);
            imagesc(IS{s,i})
%             imagesc(IS{s,i}); colormap gray; % if I want to show
%             grayscale instead of myh color transforms
            axis image
            axis off
            subplotdan(length(stains)*3,n,(s-1)*n*3 + 2*n + i);
            imagesc(err_rgb/2 + 0.5) % error between -1 and 1, put it between -0.5 and 0.5, then add 0.5 to center
            axis image
            axis off  
          %}
            %                 drawnow;
        end
    end
    
    EM = sum(EMsi(:));
    EMsave = [EMsave,EM];
    E = EM + ER;
    Esave = [Esave,E];
    if it > 1 && Esave(end) > Esave(end-1)
        
        step_mult = step_mult*2/(1 + sqrt(5));
        if step_mult < step_mult_min
            step_mult = step_mult_min;
        end
        disp(['reducing step size multiplier to ' num2str(step_mult)])
    end
    disp(['Total cost is ' num2str(E,'%g')])
    disp(['Matching cost is ' num2str(EM,'%g')])
    disp(['Regularization cost is ' num2str(ER,'%g')])
    disp(['step mult is ' num2str(step_mult)])
    % I saved these after my first run
    % save('TJ.mat','TJ')
    % save('RJ.mat','RJ')
    %{
    if it==1
        saveas(framefignum,[output_prefix 'algorithm0.png'])
    end
    if any(it==frames_to_save)
        frame = [frame,getframe(framefignum)];
        frame2Gif(frame,[output_prefix 'algorithm.gif'])    
        saveas(framefignum,[output_prefix 'algorithm.png'])
    end
    %}
    
    %%
    % plot some parameters
    danfigure(paramfignum);
    subplot(2,3,1)
    plot([Esave;EMsave;ERsave]')    
    title Matching
    
    subplot(2,3,2)
    plot(tsave')
    title Translation
    legend('x','y','z')

    % theta
    theta = acos( (trace(RM) - 1 )/2 );
    thetasave = [thetasave,theta];
    tsave = [tsave,TM];    
    subplot(2,3,3)
    plot(thetasave)
    title('Rotation angle')
    
    sxysave = [sxysave,ASxy(1,1)];
    subplot(2,3,4)
    plot(sxysave)
    title('pixel scale')
    
    szsave = [szsave,ASz(3,3)];    
    subplot(2,3,5)
    plot(szsave)
    title('slice scale')
    
    subplot(2,3,6)
    plot(ERsave)
    title('Regularization')
    
    saveas(paramfignum,[output_prefix 'params.png'])
    
    
    
    %%
    % now for scale
    % we're going to use a very similar approach
    xygrad = 0;
    deltaA = [1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,0];
    for s = 1 : length(stains)
        for i = 1 : n
            err = IS{s,i} - J{s,i};
            %tmp = reshape( DDfun(IS{s,i}(:)) * coeffs{s,i} , size(fIS{s,i}));
            %err = sum(err .* tmp, 3);
            % image gradient
            [I_x,I_y] = gradient(IS{s,i},dxJ(1),dxJ(2));
            
            
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];
            A = AJ * ASxy;
            gradmat = A*(deltaA/A);
            Xs = gradmat(1,1)*XJ{s,i} + gradmat(1,2)*YJ{s,i} + gradmat(1,3)*zJ(i) + gradmat(1,4);
            Ys = gradmat(2,1)*XJ{s,i} + gradmat(2,2)*YJ{s,i} + gradmat(2,3)*zJ(i) + gradmat(2,4);
            %         Zs = gradmat(3,1)*XJ{s,i} + gradmat(3,2)*YJ{s,i} + gradmat(3,3)*zJ(i) + gradmat(3,4);
            %         Os = gradmat(4,1)*XJ{s,i} + gradmat(4,2)*YJ{s,i} + gradmat(4,3)*zJ(i) + gradmat(4,4);
            xygrad = xygrad + -sum(sum(  err .* (I_x.*Xs + I_y.*Ys)  ))*prod(dxJ);
            
        end
    end
    
    ASxy_right_grad_mult = expm(-step_mult*epSxy * xygrad *deltaA);
    % ASxy
    disp(['xy scale is ' num2str(ASxy(1,1))])

    
    %%
    % next would be some within slice diffeo, I'm skipping this
    
    %%
    % next is the slice spacing
    zgrad = 0;
    deltaA = [0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,0];
    
    % for gradient, start after diffeo
    [phiIM_x,phiIM_y,phiIM_z] = gradient(phiIM,dxM(1),dxM(2),dxM(3));
    
    
    for s = 1 : length(stains)
        for i = 1 : n
            err = IS{s,i} - J{s,i};
            %tmp = reshape( DDfun(IS{s,i}(:)) * coeffs{s,i} , size(fIS{s,i}));
            %err = sum(err .* tmp, 3);
            % image gradient
            %         [I_x,I_y] = gradient(IS{s,i},dxJ(1),dxJ(2));
            % I can't calcluate a z gradient in this space, I have to use
            % previous space
            % let's use after the diffeo (I could have used before but that's going to get tough)
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];
            A =   AJ * ASxy * ASz * AR;
            B = inv(A);
            Xs = B(1,1)*XJ{s,i} + B(1,2)*YJ{s,i} + B(1,3)*zJ(i) + B(1,4);
            Ys = B(2,1)*XJ{s,i} + B(2,2)*YJ{s,i} + B(2,3)*zJ(i) + B(2,4);
            Zs = B(3,1)*XJ{s,i} + B(3,2)*YJ{s,i} + B(3,3)*zJ(i) + B(3,4);
            F = griddedInterpolant({yM,xM,zM},phiIM_x,'linear','nearest'); % should this be zero not nearest?
            AphiIM_x = F(Ys,Xs,Zs);
            F = griddedInterpolant({yM,xM,zM},phiIM_y,'linear','nearest'); % should this be zero not nearest?
            AphiIM_y = F(Ys,Xs,Zs);
            F = griddedInterpolant({yM,xM,zM},phiIM_z,'linear','nearest'); % should this be zero not nearest?
            AphiIM_z = F(Ys,Xs,Zs);
            
            
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];

            A = AJ * ASxy * ASz;
            gradmat = AR\(deltaA/A); % note difference from above because of different gradient! see derivation
            Xs = gradmat(1,1)*XJ{s,i} + gradmat(1,2)*YJ{s,i} + gradmat(1,3)*zJ(i) + gradmat(1,4);
            Ys = gradmat(2,1)*XJ{s,i} + gradmat(2,2)*YJ{s,i} + gradmat(2,3)*zJ(i) + gradmat(2,4);
            Zs = gradmat(3,1)*XJ{s,i} + gradmat(3,2)*YJ{s,i} + gradmat(3,3)*zJ(i) + gradmat(3,4);
%             Os = gradmat(4,1)*XJ{s,i} + gradmat(4,2)*YJ{s,i} + gradmat(4,3)*zJ(i) + gradmat(4,4); % this term is definitely zero
            zgrad = zgrad + -sum(sum(  err .* (AphiIM_x.*Xs + AphiIM_y.*Ys + AphiIM_z.*Zs)  ))*prod(dxJ);
%             zgrad = zgrad + -sum(sum(  err .* (AphiIM_z.*Zs)  ))*prod(dxJ); % I thought other terms were zero but they are not!
            
        end
    end
    
    ASz_right_grad_mult = expm(-step_mult*epSz * zgrad *deltaA);
    disp(['z scale is ' num2str(ASz(3,3))])

    
    
    %%
    % now the rigid motion of the template
    rotgrad = zeros(3,3);
    tgrad = zeros(3,1);
    for s = 1 : length(stains)
        for i = 1 : n
            err = IS{s,i} - J{s,i};
            %tmp = reshape( DDfun(IS{s,i}(:)) * coeffs{s,i} , size(fIS{s,i}));
            %err = sum(err .* tmp, 3);
            % image gradient
            %         [I_x,I_y] = gradient(IS{s,i},dxJ(1),dxJ(2));
            % I can't calcluate a z gradient in this space, I have to use
            % previous space
            % let's use after the diffeo (I could have used before but that's going to get tough)
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];
            A =   AJ * ASxy * ASz * AR;
            B = inv(A);
            Xs = B(1,1)*XJ{s,i} + B(1,2)*YJ{s,i} + B(1,3)*zJ(i) + B(1,4);
            Ys = B(2,1)*XJ{s,i} + B(2,2)*YJ{s,i} + B(2,3)*zJ(i) + B(2,4);
            Zs = B(3,1)*XJ{s,i} + B(3,2)*YJ{s,i} + B(3,3)*zJ(i) + B(3,4);
            F = griddedInterpolant({yM,xM,zM},phiIM_x,'linear','nearest'); % should this be zero not nearest?
            AphiIM_x = F(Ys,Xs,Zs);
            F = griddedInterpolant({yM,xM,zM},phiIM_y,'linear','nearest'); % should this be zero not nearest?
            AphiIM_y = F(Ys,Xs,Zs);
            F = griddedInterpolant({yM,xM,zM},phiIM_z,'linear','nearest'); % should this be zero not nearest?
            AphiIM_z = F(Ys,Xs,Zs);
            
            
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];

            A = AJ * ASxy * ASz * AR;
            
            
            rgrad = zeros(3,4);
            for i0 = 1 : 3
                for i1 = 1 : 4
                    deltaA = [double((1:4)==i0)' * double((1:4)==i1)];
                    gradmat = (deltaA/A); % again different! see derivation
                    Xs = gradmat(1,1)*XJ{s,i} + gradmat(1,2)*YJ{s,i} + gradmat(1,3)*zJ(i) + gradmat(1,4);
                    Ys = gradmat(2,1)*XJ{s,i} + gradmat(2,2)*YJ{s,i} + gradmat(2,3)*zJ(i) + gradmat(2,4);
                    Zs = gradmat(3,1)*XJ{s,i} + gradmat(3,2)*YJ{s,i} + gradmat(3,3)*zJ(i) + gradmat(3,4);
                    rgrad(i0,i1) = -sum(sum(  err .* (AphiIM_x.*Xs + AphiIM_y.*Ys + AphiIM_z.*Zs)  ))*prod(dxJ);
                end
            end
            rotgrad = rotgrad + rgrad(1:3,1:3) - rgrad(1:3,1:3)';
            tgrad = tgrad + rgrad(1:3,4);
            %         tgrad(1:2) = 0; % forget xy motions! No using a
            %         different strategy
            
            
            
        end
    end
    
% I'd like to get rid of xy here
%     tgrad(1:2) = 0; % it doesn't matter because we mix up components
%     below
% If I modeled this as two stages, one just rotation, then one just
% translation, but for now let's do it all at once
    AR_right_grad_mult =  expm([[-step_mult*epR*rotgrad,-step_mult*epT*tgrad];[0,0,0,0]]);
    disp(['TM is ' num2str(TM(1)) ' ' num2str(TM(2)) ' ' num2str(TM(3))])
    disp(['RM is '])
    disp(RM)
    %%
    % then the diffeomorphism
    % for this I will need to pull back all the data to the mri space
    % recall that scale is important when moving to mri
    lambda1 = zeros(size(IM));
    for s = 1 : length(stains)
        for i = 1 : n
            err = IS{s,i} - J{s,i};
            %tmp = reshape( DDfun(IS{s,i}(:)) * coeffs{s,i} , size(fIS{s,i}));
            %err = sum(err .* tmp, 3);
            % I need to pull back with 
            % RJ, Sxy,Sz, RM
            %AJ = [[RJ{s,i},[0;0];0,0,1],[TJ{s,i};0]; [0,0,0,1]];
            AJ = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1]];

            A = AJ*ASxy*ASz*AR;
            % sample points
            Xs = A(1,1)*XM + A(1,2)*YM + A(1,3)*ZM + A(1,4);
            Ys = A(2,1)*XM + A(2,2)*YM + A(2,3)*ZM + A(2,4);
            Zs = A(3,1)*XM + A(3,2)*YM + A(3,3)*ZM + A(3,4);
            
            zthick = 0.25;
            z_ = [zJ(i)-zthick, zJ(i), zJ(i)+zthick]; % we will zero pad
            % now I don't think the above is right
            % I think this should be thinner than ~1mm
            % this is a slicing! the slices are thin
            
            y_ = [yJ{s,i}(1)-down, yJ{s,i}, yJ{s,i}(end)+down];
            x_ = [xJ{s,i}(1)-down, xJ{s,i}, xJ{s,i}(end)+down];
            % interpolate with zero boundary conditions
            F = griddedInterpolant({y_,x_,z_},padarray(err,[1,1,1],0,'both'),'linear','nearest');
%             lambda1 = lambda1 + -F(Ys,Xs,Zs)*det(A); % should I multiply by det A?  It may be negative, how will this change things?
            % note negative sign here!
            % I'd like to just include abs(det(A))
            % block 1 has a negative jacobian, blocks 2 and 3 don't
            % block 1 is not working
            lambda1 = lambda1 + -F(Ys,Xs,Zs)*abs(det(A));
            
        end        
    end
%     figure;sliceView(lambda1)
%     return

    % I'm going to pull back the data
% I'm going to zero pad up and down, and pull back using linear
% interpolation
% this will sum to 1, and I can do this for every slice

% now we flow the error back in time
% and update v




    phi1tinvx = XM;
    phi1tinvy = YM;
    phi1tinvz = ZM;
    IMt = repmat(IM,[1,1,1,nt]);
    for t = nt : -1 : 1
        % update phi
        Xs = XM + vtx(:,:,:,t)*dt;
        Ys = YM + vty(:,:,:,t)*dt;
        Zs = ZM + vtz(:,:,:,t)*dt;
        % subtract and add identity
        F = griddedInterpolant({yM,xM,zM},phi1tinvx-XM,'linear','nearest');
        phi1tinvx = F(Ys,Xs,Zs) + Xs;
        F = griddedInterpolant({yM,xM,zM},phi1tinvy-YM,'linear','nearest');
        phi1tinvy = F(Ys,Xs,Zs) + Ys;
        F = griddedInterpolant({yM,xM,zM},phi1tinvz-ZM,'linear','nearest');
        phi1tinvz = F(Ys,Xs,Zs) + Zs;
        
        % find determinant of jacobian
        [phi1tinvx_x,phi1tinvx_y,phi1tinvx_z] = gradient(phi1tinvx,dxM(1),dxM(2),dxM(3));
        [phi1tinvy_x,phi1tinvy_y,phi1tinvy_z] = gradient(phi1tinvy,dxM(1),dxM(2),dxM(3));
        [phi1tinvz_x,phi1tinvz_y,phi1tinvz_z] = gradient(phi1tinvz,dxM(1),dxM(2),dxM(3));
        detJac = phi1tinvx_x.*(phi1tinvy_y.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_y) ...
            - phi1tinvx_y.*(phi1tinvy_x.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_x) ...
            + phi1tinvx_z.*(phi1tinvy_x.*phi1tinvz_y - phi1tinvy_y.*phi1tinvz_x);

        % update lambda
        F = griddedInterpolant({yM,xM,zM},lambda1,'linear','nearest');
        lambda = F(phi1tinvy,phi1tinvx,phi1tinvz).*detJac;
        
        % find image gradient
        [I_x,I_y,I_z] = gradient(IMt(:,:,:,t),dxM(1),dxM(2),dxM(3));
        
        % velocity field gradient gradient (note no negative sign here, it
        % is in lambda1)
        gradx = lambda.*I_x;
        grady = lambda.*I_y;
        gradz = lambda.*I_z;
        
        % smooth, and add reg
        % using vtxhat here gives opportunity for another smoothing
        % operator
        gradx = ifftn( fftn(gradx).*K + vtxhat(:,:,:,t)/sigmaR^2 , 'symmetric');
        grady = ifftn( fftn(grady).*K + vtyhat(:,:,:,t)/sigmaR^2 , 'symmetric');
        gradz = ifftn( fftn(gradz).*K + vtzhat(:,:,:,t)/sigmaR^2 , 'symmetric');
        
        
        % udpate
        vtx(:,:,:,t) = vtx(:,:,:,t) - step_mult*epV*gradx;
        vty(:,:,:,t) = vty(:,:,:,t) - step_mult*epV*grady;
        vtz(:,:,:,t) = vtz(:,:,:,t) - step_mult*epV*gradz;
        
        
    end





    
    %% 
    % do all updating now at the end
    % first the final rigids
    for s = 1 : length(stains)
        for i = 1 : n      
            % save for every 100 iterations
            % save intermediate deformations
            if mod(it,250) == 0
                %AJs = [[eye(2), [0;0];0,0,1], [TJ{s,i};0];[0,0,0,1]]; % for katie
                ARs = [RM,TM;0,0,0,1]; % R for rigid

                ASzs = diag([1,1,sz,1]);
                ASxys = diag([sxy,sxy,1,1]);

                As = ASxys * ASzs * ARs;
                Bs = inv(As);
                Bsave{i} = Bs; % not useful becasue is just for last 
            end
        end
    end
    
    
    % xy scale
    ASxy = ASxy * ASxy_right_grad_mult;
    sxy = ASxy(1,1);
    
    % z scale
    ASz = ASz * ASz_right_grad_mult;
    sz = ASz(3,3);
    
    % rigid
    AR = AR * AR_right_grad_mult;
    % this is mixing up components which I don't love
    TM = AR(1:3,4);
    RM = AR(1:3,1:3);
    
    %%    
    % save params
    save([output_prefix 'params.mat'],'RM','TM','AR','sxy','ASxy','sz','ASz','RJ','TJ')
    if ~mod(it-1,1000) || it == niter
        save([output_prefix 'deformation.mat'],'xM','yM','zM','vtx','vty','vtz','a','p')
    end
    
    
end % of iterations


% I also need a better visualization
% it would  be nice to just look at target, and do all 15 slices, and 3
% colors
% TO DO, BETTER VISUALIZATION FOR PAPER
% also contrast transforms
% I think the visualization is fine, but I need contrast transforms


%%
% for katie, remake the slices with no 2D transforms
dt = 1/nt;
phiinvx = XM;
phiinvy = YM;
phiinvz = ZM;
IMt = repmat(IM,[1,1,1,nt]);
for t = 1 : nt
    % deform image
    if t > 1
        F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
        IMt(:,:,:,t) = F(phiinvy,phiinvx,phiinvz);
    end
    % update phi
    Xs = XM - vtx(:,:,:,t)*dt;
    Ys = YM - vty(:,:,:,t)*dt;
    Zs = ZM - vtz(:,:,:,t)*dt;
    % subtract and add identity
    F = griddedInterpolant({yM,xM,zM},phiinvx-XM,'linear','nearest');
    phiinvx = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yM,xM,zM},phiinvy-YM,'linear','nearest');
    phiinvy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yM,xM,zM},phiinvz-ZM,'linear','nearest');
    phiinvz = F(Ys,Xs,Zs) + Zs;
end

F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
phiIM = F(phiinvy,phiinvx,phiinvz);


AR = [RM,TM;0,0,0,1]; % R for rigid
BR = inv(AR);

ASz = diag([1,1,sz,1]);


ASxy = diag([sxy,sxy,1,1]);



TJbar = [0;0];
thetaJbar = 0;
for i = 1 : n
    % map from R3 to R2,0 so it can be invertible
    AJ = eye(4); % no translations either 
    BJ = inv(AJ);
    danfigure(6);
        
    % I could keep pushing back
    A = AJ * ASxy * ASz * AR;
    B = inv(A);
    Bsave{i} = B;
        
    Xs = B(1,1)*XJ0_sing + B(1,2)*YJ0_sing + B(1,3)*zJ0(i) + B(1,4); % use same max and min in all slices?
    Ys = B(2,1)*XJ0_sing + B(2,2)*YJ0_sing + B(2,3)*zJ0(i) + B(2,4);
    Zs = B(3,1)*XJ0_sing + B(3,2)*YJ0_sing + B(3,3)*zJ0(i) + B(3,4);
        
    F = griddedInterpolant({yM,xM,zM},phiinvx-XM,'linear','nearest');
    phiinvBx = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yM,xM,zM},phiinvy-YM,'linear','nearest');
    phiinvBy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yM,xM,zM},phiinvz-ZM,'linear','nearest');
    phiinvBz = F(Ys,Xs,Zs) + Zs;
        
    F = griddedInterpolant({yM,xM,zM},IM,'linear','nearest');
    IS_sing{i} = F(phiinvBy,phiinvBx,phiinvBz);
    F2 = griddedInterpolant({yM,xM,zM},IMorig,'linear','nearest');
    IS_sing_orig{i} = F2(phiinvBy,phiinvBx,phiinvBz);
        
end


phisaveX = phiinvx;
phisaveY = phiinvy;
phisaveZ = phiinvz;


%%
% what do I need to save for katie?
% the MRI sliced, with no 2D transforms applied
mris = squeeze(IS_sing_orig);
save([output_prefix 'mris2.mat'],'mris')
save([output_prefix 'relTrans2.mat'],'Bsave','phisaveX', 'phisaveY','phisaveZ','sxy','sz') % have them all be single stain
end
