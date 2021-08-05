
tic;
close all; clear all;


mt_on = struct();
mt_off = struct();
fa = struct();

files = dir(fullfile(pwd, '*.nii*'));

if isempty(files)
    error('No mtOn file');
else
    
    for ii = 1:length(files)
        
        nii = files(ii).name;
        if strcmp(nii(end-1:end),'gz')
            
            if strcmp(nii(1:end-7),'mtOn')
                gunzip(nii);
                mt_on = load_untouch_nii('mtOn.nii');
            end
            
            if strcmp(nii(1:end-7),'mtOff')
                gunzip(nii);
                mt_off = load_untouch_nii('mtOff.nii');
            end
            
            
            
            
        elseif strcmp(nii(end-2:end),'nii')
            
            if strcmp(nii(1:end-4),'mtOn')
                mt_on.data = load_untouch_nii('mtOn.nii');
            end
            
            if strcmp(nii(1:end-4),'mtOff')
                mt_off.data = load_untouch_nii('mtOff.nii');
            end
            
            
        end
        
    end
end

mt_on.fprefix = mt_on.data.fileprefix;
disp(['Reading ' mt_on.fprefix ' successful.']);
temp_info=mt_on.data.hdr.hist;
disp(sprintf('Origin : %f %f %f',temp_info.qoffset_x, temp_info.qoffset_y,temp_info.qoffset_z));

temp_info=mt_on.data.hdr.hist;
mt_on.origin = zeros(1,3);
mt_on.origin(1,1) = temp_info.qoffset_x;
mt_on.origin(1,2) = temp_info.qoffset_y;
mt_on.origin(1,3) = temp_info.qoffset_z;
clear temp_info;

temp_info= mt_on.data.hdr.dime;
mt_on.spacing = zeros(1,3);
mt_on.spacing(1,1) = temp_info.pixdim(2);
mt_on.spacing(1,2) = temp_info.pixdim(3);
mt_on.spacing(1,3) = temp_info.pixdim(4);
clear temp_info;


mt_on.data.img = double(squeeze(mt_on.data.img));

mt_on.data.img = flip(mt_on.data.img,2);

mt_on.data.img = flip(mt_on.data.img,1);

mt_on.size = size(mt_on.data.img);

mt_on.direction = [1 0 0;0 1 0;0 0 1];


temp_info=mt_off.data.hdr.hist;
mt_off.origin = zeros(1,3);
mt_off.origin(1,1) = temp_info.qoffset_x;
mt_off.origin(1,2) = temp_info.qoffset_y;
mt_off.origin(1,3) = temp_info.qoffset_z;
clear temp_info;

temp_info= mt_off.data.hdr.dime;
mt_off.spacing = zeros(1,3);
mt_off.spacing(1,1) = temp_info.pixdim(2);
mt_off.spacing(1,2) = temp_info.pixdim(3);
mt_off.spacing(1,3) = temp_info.pixdim(4);
clear temp_info;

mt_off.data.img = double(squeeze(mt_off.data.img));

mt_off.data.img = flip(mt_off.data.img,2);

mt_off.data.img = flip(mt_off.data.img,1);

mt_off.size = size(mt_off.data.img);
mt_off.direction = [1 0 0;0 1 0;0 0 1];



plot = struct();

plot.slice = 25;
plot.Z = plot.slice.*ones(mt_on.size(2),mt_on.size(1));
[plot.X,plot.Y] = meshgrid(1:mt_on.size(1),1:mt_on.size(2));
plot.origin = mt_on.origin;
plot.spacing = mt_on.spacing;
plot.x = zeros(size(plot.X));
plot.y = zeros(size(plot.Y));
plot.z = zeros(size(plot.Z));
plot.direction =[1 0 0;0 1 0;0 0 1];

for i=1:mt_on.size(1)
    for j=1:mt_on.size(2)
        [outPoint]=ChangePtsCorSis(plot.origin,plot.spacing,plot.direction,[i j plot.slice],1);
        plot.x(i,j) = outPoint(1);
        plot.y(i,j) = outPoint(2);
        plot.z(i,j) = outPoint(3);
    end
end

mt_on.intensity = mat2gray(mt_on.data.img);

figure();
surf(plot.x,plot.y,plot.z,mt_on.intensity(:,:,plot.slice));
hold on;

files = dir(fullfile(pwd, '*.trk*'));

if isempty(files)
    error('No tractography input');
else
    if length(files) == 1
        trck = files.name;
        if strcmp(trck(end-2:end),'gz')
            
            gunzip(trck);
            trck = trck(1:end-4);
            [hdr,trcks] = trk_read(trck);
        else
            
            [hdr,trcks] = trk_read(trck);
            
        end
    else
        
        
        for ii=1:length(files)
            
            trck = files(ii).name;
            
            if strcmp(trck(end-2:end),'trk')
                disp(trck);
                [hdr,trcks] = trk_read(trck);
                break;
            end
            
        end
        
    end
end

disp([num2str(length(trcks)) ' tracts has been loaded']);

if isfile('transform.mat')
    transform = load('transform.mat');
    transform = transform.transform;
else
    fid = fopen('transform.txt','rt');
    transform = textscan(fid, '%f');
    fclose(fid);
    transform = reshape(cell2mat(transform),[4,4])';
end

trcks_fa = trcks;

transform(1,4) = transform(1,4)*(1/transform(1,1));
transform(2,4) = transform(2,4)*(1/transform(2,2));
transform(3,4) = transform(3,4)*(1/transform(3,3));

transform(1,1) = 1;
transform(2,2) = 1;
transform(3,3) = 1;

for ii=1:length(trcks)
    for jj = 1:length(trcks(ii).matrix)
        trcks(ii).matrix(jj,:) = trcks(ii).matrix(jj,:)*transform(1:3,1:3) + mt_on.origin - transform(1:3,4)';
    end
end

trk_plot(hdr,trcks)

trcks_fa = trcks;

for ii=1:length(trcks)
    for jj = 1:length(trcks(ii).matrix)
        trcks(ii).matrix(jj,:) = trcks(ii).matrix(jj,:)*transform(1:3,1:3) + mt_on.origin - transform(1:3,4)';
    end
end


MTR = mt_on;
MTR.data.img = abs(mt_on.data.img - mt_off.data.img)./mt_off.data.img;
MTR.data.img(isinf(MTR.data.img)) = 0;
MTR.data.img(isnan(MTR.data.img)) = 0;
MTR.data.img(MTR.data.img>1) = 1;

disp('Interpolating mtr');
[IJKneighbors,interpWeighs] = idw3dInterp(trcks,MTR);

for ii=1:length(trcks)
    for jj=1:length(trcks(ii).matrix)

        idx = squeeze(IJKneighbors(ii).mat(jj,:,:));

tmp = [];
for kk=1:8
    
    tmp = [tmp MTR.data.img(idx(1,kk),idx(2,kk),idx(3,kk))];
end

trcks(ii).matrix(jj,4) = sum(interpWeighs(ii).mat(jj,:).*tmp);



    end
end

figure();
hdr.n_scalars = 1;
surf(plot.x,plot.y,plot.z,mt_on.intensity(:,:,plot.slice));
hold on;
trk_plot(hdr,trcks,[],[],'scalar',1)
toc;


