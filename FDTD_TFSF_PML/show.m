clc;clear all;close all;
disp('正在计算，请稍等...');
system('FDTD3D_TFSF.exe');

frame = 0;
basename1 = './data/plane-z';
filename1 = sprintf('%s.%d', basename1, frame);
fid1 = fopen(filename1, 'rb');
basename2 = './data/plane-y';
filename2 = sprintf('%s.%d', basename2, frame);
fid2 = fopen(filename2, 'rb');
if fid1 == -1
    error(['raw2movie: initial frame not found: ', filename1])
end
if fid2 == -1
    error(['raw2movie: initial frame not found: ', filename2])
end

while fid1 ~= -1
    timestep=int2str(frame + 1);
    size_x1 = fread(fid1, 1, 'single');
    size_y1 = fread(fid1, 1, 'single');
    dataz = flipud(transpose(reshape(fread(fid1, size_x1 * size_y1, 'single'), size_x1, size_y1)));
    size_x2 = fread(fid2, 1, 'single');
    size_z2 = fread(fid2, 1, 'single');
    datay = flipud(transpose(reshape(fread(fid2, size_x2 * size_z2, 'single'), size_x2, size_z2)));

    subplot(211),imagesc(dataz');
    shading flat;
    caxis([-0.5 0.5]); 
    colorbar;
    axis image; axis xy; 
    title(['Ez(x,y,z=50), time step = ',timestep]);
    xlabel('x coordinate'); ylabel('y coordinate');
    
    subplot(212),imagesc(datay');
    shading flat;
    caxis([-0.5 0.5]); 
    colorbar;
    axis image; axis xy;
    title(['Ez(x,y=50,z), time step = ',timestep]);
    xlabel('x coordinate'); ylabel('z coordinate');

    pause(0.001)
    fclose(fid1);
    fclose(fid2);
    frame = frame + 1;
    filename1 = sprintf('%s.%d', basename1, frame);
    fid1 = fopen(filename1, 'rb');
    filename2 = sprintf('%s.%d', basename2, frame);
    fid2 = fopen(filename2, 'rb');
end