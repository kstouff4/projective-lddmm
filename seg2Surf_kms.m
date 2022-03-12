%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use delaunay triangulation to create surfaces from segmentations
%
%  Input should be a bunch of filenames, the last arguments are the minimum
%  edge size, resolution, origin as -1 = corner and 1 = center
%
%  why do it this way?  Because bash is terrible and I can't figure out how
%  to pass an asterisk
%  
%

function [F,V] = seg2Surf_kms(varargin)
disp('inputs')
disp(varargin)
% keyboard
% test
% pattern = '/cis/home/dtward/Documents/biocard/filteringStudy/targetData_v02/0186193_1_6_hippo_*_daniel1.img';



% files = dir(pattern);
% [d,f,e] = fileparts(pattern);
% filenames = {};
% disp(['Files that match pattern ' pattern ':'])
% for i = 1 : length(files)
%     filenames{i} = [d filesep files(i).name];
%     disp(filenames{i})
% end
filenames = varargin;
if ~exist(filenames{end-2},'file')
    minD = str2num(filenames{end-2});
    res = str2num(filenames{end-1});
    orig = str2num(filenames{end});
    filenames = filenames(1:end-3);
else
    minD = NaN;
    res = str2num(filenames{end-1});
    orig = str2num(filenames{end});
    filenames = filenames(1:end-2);
end


draw = 0;

% add path for reading and writing analyze images
addpath /cis/home/dtward/Functions/avwQuiet
% add path for reading and writing byu surfaces
addpath /cis/home/dtward/Functions/byu
% add path for restricted delaunay triangulation
addpath /cis/home/dtward/Documents/sandbox/delaunay

% loop through filenames
for i = 1 : length(filenames)
    disp(['Starting ' filenames{i}])
    if draw
        clf;
    end
    % check for img extension
    if ~strcmp( filenames{i}(end-3:end) , '.img')
        warning([filenames{i} ' is not analyze img, skipping'])
        continue
    end
    if ~exist( filenames{i}, 'file')
        warning([filenames{i} ' does not exist, skipping']);
        continue
    end
    
    % load it
    avw = avw_img_read(filenames{i},0);
    
    
    % extract labels
    labels = unique(avw.img);
    for l = 1 : length(labels)
    
        if labels(l) == 0
            % 0 is always background, we don't want this
            continue
        end
        
        % if there are many labels, this is probably a smooth image
        segMax = 100;
        if length(labels) > segMax
            warning(['More than ' num2str(segMax) ' labels, interpretting as non-binary segmentation.'])
            I = avw.img;
        else % otherwise pick this label
            I = double( avw.img == labels(l) );
        end
        
        % get the size
        nx = double(avw.hdr.dime.dim([3,2,4]));
        dx = double(avw.hdr.dime.pixdim([3,2,4])); % ignore pixel dimension here
        disp('Before and after dx')
        disp(dx)
        if (res > 0)
            disp('changing res')
            dx = dx./dx;
            dx = dx.*res; % assume 0.125 mm
        end
        disp(dx)
        meanx = mean((0:nx(1)-1)*dx(1))
        meany = mean((0:nx(2)-1)*dx(2))
        meanz = mean((0:nx(3)-1)*dx(3))

        
        % first thing is subcube, with 3 voxel buffer
        zsum = squeeze( sum(sum(I,1),2) );
        zmin = find(zsum, 1, 'first');
        zmax = find(zsum, 1, 'last');
        ysum = squeeze( sum(sum(I,2),3) );
        ymin = find(ysum,1,'first');
        ymax = find(ysum,1,'last');
        xsum = squeeze( sum(sum(I,1),3));
        xmin = find(xsum,1,'first');
        xmax = find(xsum,1,'last');
        
        xmin = max(xmin-3, 1);
        ymin = max(ymin-3, 1);
        zmin = max(zmin-3, 1);
        
        xmax = min(xmax+3,nx(1));
        ymax = min(ymax+3,nx(2));
        zmax = min(zmax+3,nx(3));
        
        I = I(ymin:ymax,xmin:xmax,zmin:zmax);
        x0 = ([xmin,ymin,zmin]-1) .* dx; % if xmin, the index, is 1, we want x0 to be 0
        nx = [size(I,2),size(I,1),size(I,3)];
        x = (0:nx(1)-1)*dx(1) + x0(1);
        y = (0:nx(2)-1)*dx(2) + x0(2);
        z = (0:nx(3)-1)*dx(3) + x0(3);
        
        % smooth it by 1 voxel
        disp('About to smooth')
        smooth = 1;
%         % try less smoothing
%         smooth = 0.5;
        x_ = -3 : 3;
        y_  = -3 : 3;        
        z_ = -3 : 3;
        [X_,Y_,Z_] = meshgrid(x_,y_,z_);
        K = exp(-(X_.^2 + Y_.^2 + Z_.^2)/2/smooth^2);
        K = K/sum(K(:));
        I = convn(I,K,'same');
        
        % upsample?
        % let's go to uniform, twice resolution of smallest dx        
        upsample = 1;
        if upsample
            disp(['about to upsample'])
            lx = nx.*dx;
            dx_ = min(dx)*[1 1 1]/2;
            
            x_ = (0:dx_(1):lx(1)) + x0(1);
            y_ = (0:dx_(2):lx(2)) + x0(2);
            z_ = (0:dx_(3):lx(3)) + x0(3);
            [X_,Y_,Z_] = meshgrid(x_,y_,z_);
            [X,Y,Z] = meshgrid(x,y,z);
            I = interp3(X,Y,Z,I,X_,Y_,Z_,'linear',0);
            dx = dx_;
            nx = [size(I,2),size(I,1),size(I,3)];
            x = x_;
            y = y_;
            z = z_;
        end
        
%         keyboard
        % triangulate
        disp('about to create isosurface')
        OPT = struct;
        OPT.dx = dx;
        [f,v] = isosurface(x,y,z,I,max(I(:))/2);
        if draw
            patch('faces',f,'vertices',v,'facecolor','r','edgecolor','k','facealpha',0.5);
        end
        OPT.nIdeal = size(v,1);  % I will definitely have less vertices than this
%         OPT.nIdeal = min(OPT.nIdeal,5000); % a hack
        OPT.minD = sqrt(sum(dx.^2));
        if ~isnan(minD)
            if ischar(minD)
                OPT.minD = str2num( minD ); % use it if specified
            else
                OPT.minD = minD;
            end
        else
            warning(['Minimum edge length not specified, using ' num2str(OPT.minD)]);
        end
        disp('about to triangulate')
%         OPT.verbose = 1;
        [F,V] = restrictedDelaunayFromImage(I,OPT);
        V = bsxfun(@plus, V, x0);
        if draw
            patch('faces',F,'vertices',V,'facecolor','b','edgecolor','k','facealpha',0.5);
            view(30,30)
            axis image
            light
        end
        
        if (orig > 0)
            disp('Means of shift are')
            disp(meanx)
            disp(meany)
            disp(meanz)
            V(:,1) = V(:,1) - meanx;
            V(:,2) = V(:,2) - meany;
            V(:,3) = V(:,3) - meanz;
        end
        
        % write out, locally, with the same filename, and extension added
        disp('about to write out')
        [d,f,e] = fileparts(filenames{i});
        disp('file path is ' )
        disp(d)
        disp(f)
        
        if length(labels) > segMax
            disp(['Writing file ' f '.byu'])
            writebyu(size(V,1), size(F,1), size(F,1)*3, V, F, [f '.byu'])
            break; % only once
        else
            disp(['Writing file ' f '_label' num2str(labels(l),'%02d') '.byu'])
            writebyu(size(V,1), size(F,1), size(F,1)*3, V, F, [f '_label' num2str(labels(l),'%02d') '.byu']);
        end
    end
    
end
