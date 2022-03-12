function session2Analyze(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% session2Analyze(filename)
%
%    filename is the name of a Seg3D session.  either .ses (version 1) or
%    .s3d (version2)
%
%    The program will output one analyze image for each layer in the
%    current directory.  This includes grayscale images and segmentations.  
%    The name will be the name of the layer.
%
%    TO DO: If layer names are not unique they will be overwrited
%
%    Depenencies: analyze reading and writing in dtward's function
%    directory  /cis/home/dtward/Functions/avwQuiet
%
%    Permissions:  This program will need to write to the current
%    directory, AND generally to the /tmp directory.
%

% test examples
% filename = '/cis/project/adni/work/timothy/eliza_41_non_bifurcated_less_10_high_atrophy/CIS_ADNI/CS_QC_014_S_4263_moBL.seg3dproj/CS_QC_014_S_4263_moBL.s3d'
% filename = '/cis/project/biocard/work/chelsea/CS_0224447_1_4_up2TL.ses';


% add analyze directory
addpath /cis/home/dtward/Functions/avwQuiet

% parse input
if ~exist(filename,'file')
    error(['Input filename ' filename ' does not exist']);
end
disp(['Processing session file ' filename])
[directory,file,extension] = fileparts(filename);

% check extension
version = 0;
if strcmp(extension,'.s3d')
    disp('Found Seg3D 2 session')
    version = 2;
elseif strcmp(extension,'.ses')
    disp('Found Seg3D 1 session')
    version = 1;
end
if version == 0
    error('Cannot recognize Seg3D version')
end


if version == 2
    % first thing is to load xml file
    sessionsDir = [directory filesep 'sessions' filesep];
    if ~exist(sessionsDir,'dir')
        error('Session subdirectory does not exist');
    end
    xmlFiles = dir([sessionsDir '*.xml']);
    if isempty(xmlFiles)
        error('No session xml file found')
    end
    if length(xmlFiles) > 1
%         error('More than one session xml file not supported')
        xmlFiles.name
        warning('More than one session xml file not supported, trying last file by date')
        % or try this
        % need last by date
%         xmlFiles = xmlFiles(end);
        dates = cell2mat({xmlFiles.datenum});
        xmlFiles = xmlFiles(dates == max(dates));        
    end
    
    
    
    
    doc = xmlread([sessionsDir xmlFiles(1).name]);
    % I can call fields(doc) and methods(doc)
    layers = doc.getElementsByTagName('layers').item(0);
    % this returns some kind of list, get the first element
    % the items in layers are now its children    
%     layer = layers.getChildNodes(); % this is redundant, layer == layers
    % loop through
    
    layersNumber = {};
    layersType = {};
    layersName = {};
    layersGeneration = {};
    layersBit = {};
    for i = 1 : layers.getLength()
        layerNumber = '-';
        layerType = '-';
        layerGeneration = '-';
        layerBit = '-';
        layerName = '-';
        
        node = layers.item(i-1);
        if node.getNodeType == doc.TEXT_NODE
            % there seems to be a bunch of dashes in there for some reason
            continue
        end
        % otherwise let's get the layer namne
        layerNumber = char( node.getNodeName );
        % get the type
        layerType = char( node.getAttribute('type') );
        % more data is in the children
        states = node.getChildNodes();
        for j = 1 : states.getLength()
            state = states.item(j-1);
            if state.getNodeType == doc.TEXT_NODE
                % there seems to be a bunch of dashes in there for some reason
                continue
            end
            % we want generation, this tells us the filename to read
            attribute = char( state.getAttribute('id') );
            data = char( state.item(0).getData );
            if strcmp( attribute, 'generation' )                
                layerGeneration = data; % why use item? I don't know
            end
            
            % name will tell us what to call the output
            if strcmp( attribute, 'name' )
                layerName = data;
            end
            
            % I will also need bit
            if strcmp(attribute,'bit')
                layerBit = data;
            end
            
            
        end
        
        layersNumber = [layersNumber,layerNumber];
        layersType = [layersType,layerType];
        layersName = [layersName,layerName];
        layersGeneration = [layersGeneration,layerGeneration];
        layersBit = [layersBit,layerBit];
        
        
    end
    
    
    
    
    
    % now we'll loop through all the layers and extract data one by one
    for i = 1 : length(layersNumber)
        % find directory
        dataDir = [directory filesep 'data' filesep];
        if ~exist(dataDir,'dir')
            error('Data subdirectory does not exist');
        end
        
        % find data
        nrrdFile = [dataDir layersGeneration{i} '.nrrd'];
        if ~exist(nrrdFile,'file')
            error(['NRRD file ' nrrdFile ' does not exist']);
        end
        
        % open the file
        [fid,err] = fopen(nrrdFile,'rb');
        if fid == -1
            error(['Problem opening nrrd file ' nrrdFile ': ' err]);
        end
        
        % read lines of text until you find a blank line
        zipped = 0;
        dx = [0.0 0.0 0.0];
        nx = [0 0 0];
        type = '';
        while 1
            line = fgetl(fid);
            if isempty(line)
                break;
            end
            % get data I need
            [T,R] = strtok(line,':');
            if length(R) > 1
                value = strtrim(R(2:end));
            end
            if strfind(line,'type')
                type = value;
            end
            if strfind(line,'sizes');
                nx = sscanf(value,'%d %d %d')';
            end
            if strfind(line,'spacings');
                dx = sscanf(value,'%f %f %f')';
            end
            if strfind(line,'encoding')
                if strfind(line,'gzip')
                    zipped = 1;                
                end
            end
        end
        
        % now get the data
        if zipped
            data = fread(fid);
            fclose(fid);
            % write it to tmp
            
            fid = fopen('/tmp/tmp.bin','wb');
            fwrite(fid,data);
            fclose(fid);
            % unzip it
            gunzip('/tmp/tmp.bin');
            fid = fopen('/tmp/tmp','rb');
        end
        
        % set up analyze
        avw = avw_hdr_make();
        avw.hdr.dime.dim([3,2,4]) = nx;
        avw.hdr.dime.pixdim([3,2,4]) = dx;
        % /*Acceptable values for datatype are*/
        % #define DT_NONE             0
        % #define DT_UNKNOWN          0    /*Unknown data type*/
        % #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/
        % #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/
        % #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/
        % #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/
        % #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/
        % #define DT_COMPLEX         32    /*Complex (64 bits per voxel; 2 floating point numbers)/*
        % #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/
        % #define DT_RGB            128    /*A Red-Green-Blue datatype*/
        % #define DT_ALL            255    /*Undocumented*/
        
        if strcmp(type,'short')
            datatype = 'int16';
            avw.hdr.dime.bitpix = 16;
            avw.hdr.dime.datatype = 4;
        elseif strcmp(type,'unsigned char')
            datatype = 'uint8';
            avw.hdr.dime.bitpix = 8;
            avw.hdr.dime.datatype = 2;
        elseif strcmp(type,'float')
            datatype = 'float32';
            avw.hdr.dime.bitpix = 32;
            avw.hdr.dime.datatype = 16;
        end
        data = fread(fid,datatype);
        fclose(fid);
        data = permute( reshape(data,nx([1,2,3])) , [2,1,3]); 
        if strcmp(layersType{i},'mask')
            data = double(bitget(uint8(data),str2num(layersBit{i})+1));
        end
        % I don't know if this is correct, I'll have to check with an image where the first two sizes are not the same
        % yes this is correct with 1,2,3
        avw.img = data;
%         avw_img_write(avw,[layersNumber{i} '_' layersName{i} '.img']);
        avw_img_write(avw,[layersName{i} '.img'],0);
        % don't do layer number for consistency with below
        
    end
    
elseif version == 1
    % uncompress this file to tmp
    randstr = ['session2Analyze-'];
    for i = 1 : 32% I had done 16 but I was getting collisions!
        randstr = [randstr,num2str(randi(100)-1)];
    end
    randstr = [randstr '/'];
    
    tmpdir = ['/tmp/' randstr];
    mkdir(tmpdir);
    gunzip(filename , tmpdir);
    myfile = dir(tmpdir);
    untar([tmpdir myfile(3).name],tmpdir);
    % now there should be an xml file
    xmlfile = dir([tmpdir '*.xml']);
    if isempty(xmlfile)
        error('No xml file found in session')
    end
    if length(xmlfile) > 1
        error('Multiple xml files not supported');
    end
    doc = xmlread([tmpdir xmlfile.name]);
    volumes = doc.getElementsByTagName('volume');
    layersNumber = {};
    layersType = {};
    layersName = {};
    layersGeneration = {};
    layersBit = {};
    
    for i = 1 : volumes.getLength()
        layerNumber = '-';
        layerType = '-';
        layerGeneration = '-';
        layerBit = '-';
        layerName = '-';
        
        volume = volumes.item(i-1);
        for j = 1 : volume.getLength()
            var = volume.item(j-1);
            if var.getNodeType() == doc.TEXT_NODE
                % some dashes and things that are not necessary
%                 disp('text node data:')
%                 disp(var.getData)
                continue
            end
            
            attribute = char( var.getAttribute('name') );
            
%             var.getLength()
            try
                value = char( var.item(0).getData() );
            catch
                % I should just skip values that I can't get
                continue
            end
            
            if strcmp(attribute,'name')
                layerName = value;
            end
            if strcmp(attribute,'filename')
                layerGeneration = value;
                if ~isempty(strfind(layerGeneration,'data'))
                    layerType = 'data';
                else
                    layerType = 'mask';
                end
            end
            if strcmp(attribute,'label')
                layerBit = value;
            end
            
            
        end
        
        layersNumber = [layersNumber,layerNumber];
        layersType = [layersType,layerType];
        layersName = [layersName,layerName];
        layersGeneration = [layersGeneration,layerGeneration];
        layersBit = [layersBit,layerBit];
        
    end
    
        
%     keyboard
    
    % now we'll loop through all the volumes and extract data one by one
    for i = 1 : length(layersNumber)        
        % find data
        nrrdFile = [tmpdir layersGeneration{i}];
        if ~exist(nrrdFile,'file')
            error(['NRRD file ' nrrdFile ' does not exist']);
        end
        
        % open the file
        [fid,err] = fopen(nrrdFile,'rb');
        if fid == -1
            error(['Problem opening nrrd file ' nrrdFile ': ' err]);
        end
        
        % read lines of text until you find a blank line
        zipped = 0;
        dx = [0.0 0.0 0.0];
        nx = [0 0 0];
        type = '';
        while 1
            line = fgetl(fid);
            if isempty(line)
                break;
            end
            % get data I need
            [T,R] = strtok(line,':');
            if length(R) > 1
                value = strtrim(R(2:end));
            end
            if strfind(line,'type')
                type = value;
            end
            if strfind(line,'sizes');
                nx = sscanf(value,'%d %d %d')';
            end
            if strfind(line,'spacings');
                dx = sscanf(value,'%f %f %f')';
            end
            if strfind(line,'encoding')
                if strfind(line,'gzip')
                    zipped = 1;                
                end
            end
        end
        
        % now get the data
        if zipped
            data = fread(fid);
            fclose(fid);
            % write it to tmp
            
            fid = fopen('/tmp/tmp.bin','wb');
            fwrite(fid,data);
            fclose(fid);
            % unzip it
            gunzip('/tmp/tmp.bin');
            fid = fopen('/tmp/tmp','rb');
        end
        
        % set up analyze
        avw = avw_hdr_make();
        avw.hdr.dime.dim([3,2,4]) = nx;
        avw.hdr.dime.pixdim([3,2,4]) = dx;
        % /*Acceptable values for datatype are*/
        % #define DT_NONE             0
        % #define DT_UNKNOWN          0    /*Unknown data type*/
        % #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/
        % #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/
        % #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/
        % #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/
        % #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/
        % #define DT_COMPLEX         32    /*Complex (64 bits per voxel; 2 floating point numbers)/*
        % #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/
        % #define DT_RGB            128    /*A Red-Green-Blue datatype*/
        % #define DT_ALL            255    /*Undocumented*/
        
        if strcmp(type,'short')
            datatype = 'int16';
            avw.hdr.dime.bitpix = 16;
            avw.hdr.dime.datatype = 4;
        elseif strcmp(type,'unsigned char')
            datatype = 'uint8';
            avw.hdr.dime.bitpix = 8;
            avw.hdr.dime.datatype = 2;
        elseif strcmp(type,'float')
            datatype = 'float32';
            avw.hdr.dime.bitpix = 32;
            avw.hdr.dime.datatype = 16;
        end
        data = fread(fid,datatype);
        fclose(fid);
        data = permute( reshape(data,nx([1,2,3])) , [2,1,3]); 
        if strcmp(layersType{i},'mask')
            data = double(bitget(uint8(data),log(str2num(layersBit{i}))/log(2)+1));
        end
        % I don't know if this is correct, I'll have to check with an image where the first two sizes are not the same
        % yes this is correct with 1,2,3
        avw.img = data;
%         avw_img_write(avw,[layersNumber{i} '_' layersName{i} '.img']);

        % commenting out, fix permutations
%         avw_img_write(avw,[layersName{i} '.img'],0); % no layer number here
        
        
        % note about permutations, this is NOT coming up the same way as it
        % was saved
        % fixing this on november 15, 2016
        % this will load into seg 3D the same way as the original nii image
        avw.hdr.dime.dim([2,3,4]) = nx;
        avw.hdr.dime.pixdim([2,3,4]) = dx;
        avw.img = permute(data,[2,1,3]);
        avw_img_write(avw,[layersName{i} '.img'],0); % no layer number here
        
        
    end
    rmdir(tmpdir, 's')   
end
