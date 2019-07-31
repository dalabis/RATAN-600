function ScanImageOverlay(filenameRatan, filenameImage)
%%
% RATAN-600
scan = fitsread(filenameRatan);
scan = scan(:,:,5:end-1);
infoRatan = fitsinfo(filenameRatan);
pixRatanSuccess = 0;
solarCenterRatanSuccess = 0;
solarRadiusArcsecRatanSuccess = 0;
solarRASuccess = 0;
solarDecSuccess = 0;
i = 1;

while (i <= size(infoRatan.PrimaryData.Keywords, 1)) && ~(pixRatanSuccess && solarCenterRatanSuccess && solarRadiusArcsecRatanSuccess && solarRASuccess && solarDecSuccess)  
    if isequal(infoRatan.PrimaryData.Keywords{i,1}, 'CDELT1')
        pixRatan = infoRatan.PrimaryData.Keywords{i,2};
        pixRatanSuccess = 1;
    end
    
    if isequal(infoRatan.PrimaryData.Keywords{i,1}, 'CRPIX1')
        solarCenterRatan = infoRatan.PrimaryData.Keywords{i,2};
        solarCenterRatanSuccess = 1;
    end
    
    if isequal(infoRatan.PrimaryData.Keywords{i,1}, 'SOLAR_R')
        solarRadiusArcsecRatan = infoRatan.PrimaryData.Keywords{i,2};
        solarRadiusArcsecRatanSuccess = 1;
    end
    
    if isequal(infoRatan.PrimaryData.Keywords{i,1}, 'SOL_RA')
        solarRA = infoRatan.PrimaryData.Keywords{i,2};
        solarRASuccess = 1;
    end
    
    if isequal(infoRatan.PrimaryData.Keywords{i,1}, 'SOL_DEC')
        solarDec = infoRatan.PrimaryData.Keywords{i,2};
        solarDecSuccess = 1;
    end
    
    i = i + 1;

end

%%
% Another telescope 2D image
% CORONAS SPIRIT SCT 171A
% SOHO 171A
%rotateSpirit = 104.5;
%solarRadiusPixSpirit = 393.3;
image = fitsread(filenameImage);
infoImage = fitsinfo(filenameImage);
rotateSpiritSuccess = 0;
solarRadiusPixSpiritSuccess = 0;
nameTelescopeSuccess = 0;
CRPIX1_ImageSuccess = 0;
CRPIX2_ImageSuccess = 0;
i = 1;

while (i <= size(infoImage.PrimaryData.Keywords, 1)) && ~(rotateSpiritSuccess && solarRadiusPixSpiritSuccess && nameTelescopeSuccess && CRPIX1_ImageSuccess && CRPIX2_ImageSuccess)
    if isequal(infoImage.PrimaryData.Keywords{i,1}, 'SC_ROLL')
        rotateSpirit = infoImage.PrimaryData.Keywords{i,2};
        rotateSpiritSuccess = 1;
    end
    
    if isequal(infoImage.PrimaryData.Keywords{i,1}, 'SOLAR_R')
        solarRadiusPixSpirit = infoImage.PrimaryData.Keywords{i,2};
        solarRadiusPixSpiritSuccess = 1;
    end
    
    if isequal(infoImage.PrimaryData.Keywords{i,1}, 'TELESCOP')
        nameTelescope = infoImage.PrimaryData.Keywords{i,2};
        nameTelescopeSuccess = 1;
    end
    
    if isequal(infoImage.PrimaryData.Keywords{i,1}, 'CRPIX1')
        CRPIX1_Image = infoImage.PrimaryData.Keywords{i,2};
        CRPIX1_ImageSuccess = 1;
    end
    
    if isequal(infoImage.PrimaryData.Keywords{i,1}, 'CRPIX2')
        CRPIX2_Image = infoImage.PrimaryData.Keywords{i,2};
        CRPIX2_ImageSuccess = 1;
    end
    
    i = i + 1;
    
end

%%
if isequal(nameTelescope,'SPIRIT CORONAS-F')
    figure
    imshow(image,[250 1500])
    hold on

    N = 5;
    [x,y] = ginput(N);

    % lsqr optimization
    a1 = 2/N*sum(x)^2 - 2*sum(x.^2);
    b1 = 2/N*sum(x)*sum(y) - 2*sum(x.*y);
    %a2 = b1;
    b2 = 2/N*sum(y)^2 - 2*sum(y.^2);
    c1 = sum(x.^3) + sum(x.*y.^2) - 1/N*sum(x)*sum(x.^2+y.^2);
    c2 = sum(y.^3) + sum(y.*x.^2) - 1/N*sum(y)*sum(x.^2+y.^2);

    x0 = (b1*c2 - c1*b2) / (a1*b2 - b1^2);
    y0 = (b1*c1 - a1*c2) / (a1*b2 - b1^2);

    r = solarRadiusPixSpirit;
    c = [x0 y0];
    pos = [c-r 2*r 2*r];
    rectangle('Position',pos,'Curvature',[1 1])
else
    x0 = CRPIX1_Image;
    y0 = CRPIX2_Image;
end

%%
newIm(1:fix(4*solarRadiusPixSpirit),1:fix(4*solarRadiusPixSpirit)) = 1500;
newIm(fix(size(newIm,2)/2 - y0)+1:fix(size(newIm,2)/2 - y0)+size(image,2),...
    fix(size(newIm,1)/2 - x0)+1:fix(size(newIm,1)/2 - x0)+size(image,1)) = image;

% rotate image
if isequal(nameTelescope,'SPIRIT CORONAS-F')
    newIm = imrotate(newIm, (180 + rotateSpirit - atan(sin(solarDec/180*pi) * tan((solarRA-90)/180*pi))/pi*180),'crop');
else
    newIm = imrotate(newIm, (rotateSpirit + atan(sin(solarDec/180*pi) * tan((solarRA-90)/180*pi))/pi*180), 'crop');
    newIm = flipud(newIm);
end

%%
figure('units','centimeters','position',[2 2 9 9])
subplot('position', [1/9 1/9 7/9 7/9])
imagesc(newIm,[250 1500])
colormap(gray)
axis on
hold on

xlim([ (size(newIm,2)/2 - 1.3*solarRadiusPixSpirit)  (size(newIm,2)/2 + 1.3*solarRadiusPixSpirit) ])
ylim([ (size(newIm,1)/2 - 1.3*solarRadiusPixSpirit)  (size(newIm,1)/2 + 1.3*solarRadiusPixSpirit) ])

%%%
% unknown parameter
unknownSwing = -20;
%%%

xi = linspace(1, size(scan,2), size(scan,2));
xi = xi - (solarCenterRatan + unknownSwing - solarRadiusArcsecRatan / pixRatan);
xi = xi ./ (solarRadiusArcsecRatan / pixRatan) .* solarRadiusPixSpirit;
xi = xi + (size(newIm,1)/2 - solarRadiusPixSpirit);

minScan = min(min(scan(1,:,:) - scan(2,:,:)));
maxScan = max(max(scan(1,:,:) + scan(2,:,:)));
deltaPix = 2*1.3*solarRadiusPixSpirit / (maxScan - minScan);
swingPix = 3.3*solarRadiusPixSpirit + minScan * deltaPix;

scanDraw = scan(:,:,3:end);
scanDraw(1,:,:) = -((scan(1,:,3:end) + scan(2,:,3:end))).*deltaPix + swingPix;
scanDraw(2,:,:) = -((scan(1,:,3:end) - scan(2,:,3:end))).*deltaPix + swingPix;

color = jet(size(scanDraw,3));

for i = 1:size(scanDraw,3)
    plot(xi,scanDraw(1,:,i),'color',color(i,:))
    plot(xi,scanDraw(2,:,i),'color',color(i,:))
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'xtick', [-1200;-800;-400;0;400;800;1200] .* solarRadiusPixSpirit ./ (solarRadiusArcsecRatan / pixRatan) ./ pixRatan + size(newIm,2)/2)
set(gca, 'xticklabel', {'-1200';'-800';'-400';'0';'400';'800';'1200'})
set(gca, 'ytick', -[6;5;4;3;2;1;0;-1].*10^4 .* deltaPix + swingPix)
set(gca, 'yticklabel', {'6';'5';'4';'3';'2';'1';'0';'-1'})

xlabel('{Distance from the solar center, arcsec}', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('{$T_a$, $10^4$K}', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 10)

line([ (size(newIm,2)/2 - solarRadiusArcsecRatan * solarRadiusPixSpirit / (solarRadiusArcsecRatan / pixRatan) / pixRatan) (size(newIm,2)/2 - solarRadiusArcsecRatan * solarRadiusPixSpirit / (solarRadiusArcsecRatan / pixRatan) / pixRatan) ],...
    [ 1 size(newIm,1) ],'color','k','linestyle','--')
line([ size(newIm,2)/2 size(newIm,2)/2 ],...
    [ 1 size(newIm,1) ],'color','k','linestyle','--')
line([ (size(newIm,2)/2 + solarRadiusArcsecRatan * solarRadiusPixSpirit / (solarRadiusArcsecRatan / pixRatan) / pixRatan) (size(newIm,2)/2 + solarRadiusArcsecRatan * solarRadiusPixSpirit / (solarRadiusArcsecRatan / pixRatan) / pixRatan) ],...
    [ 1 size(newIm,1) ],'color','k','linestyle','--')