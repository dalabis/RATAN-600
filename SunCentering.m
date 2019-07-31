clear
close all

%% read .fits file and header
fileName = '011022sun0_out_edit.fits';
% data
data = fitsread(fileName);
% deliting bad frequencies
badFreq = 2;
data = data(:,:,3:end);
% header
info = fitsinfo(fileName);
% read from header frequencies
endSuccess = 0;
commentSuccess = 0;
freqNum = 1;
counter = 1;

while ~endSuccess
    if isequal(info.PrimaryData.Keywords{counter, 1}, 'COMMENT')
        commentSuccess = 1;
    end
    
    if isequal(info.PrimaryData.Keywords{counter+1, 1}, 'END')
        endSuccess = 1;
    end
    
    counter = counter + 1;
    
    if commentSuccess == 1 && endSuccess == 0
        freq(freqNum) = info.PrimaryData.Keywords{counter, 2};
        freqNum = freqNum + 1;
    end
end

% deleting bad frequencies
freq = freq(1+badFreq:end);

% transition from pixels to arcsec units
CDELT = 4.9405;
% solar radius, pix
R = 964.2245912664 / CDELT;
% solar center, pix
CRPIX = 535.505615;
% frequencies and antenna patern values from "Work Scan" program
f = [  0.985;   1.015;   1.045;   1.670;   1.760;   1.860;   1.950; ...
       2.050;   2.150;   2.270;   2.610;   2.720;   2.830;   2.950; ...
       3.080;   3.210;   3.350;   3.670;   3.950;   4.270;   4.600;   4.950;   5.700; ...
       6.080;   6.500;   6.950;   7.350;   7.830;   8.400;   8.750;   9.350;   9.800; ...
      10.350;  10.950;  11.250;  12.950;  13.400;  14.250;  14.750;  15.650; ...
      16.400                                                      ];
D = [237.410; 235.200; 233.000; 188.000; 181.810; 174.940; 168.760; ...
     162.240; 156.070; 148.660; 128.680; 122.880; 117.090; 110.770; ...
     104.660;  99.010;  92.930;  80.500;  70.770;  61.770;  53.580;  46.550;  37.200; ...
      33.250;  31.330;  29.270;  28.560;  27.910;  27.420;  27.190;  26.730;  26.360; ...
      25.740;  24.970;  24.500;  21.580;  20.790;  19.330;  18.510;  17.290; ...
      16.360                                                      ];

%%
% founding antenna pattern values from frequencies 
% horizontal antenna pattern
Dfreq(1:length(freq)) = 0;
j = 1;

for i = 1:length(f)
    if j <= length(freq) && freq(j) == f(i)
        Dfreq(j) = D(i);
        j = j + 1;
    end
end

% vertical antenna pattern
DVertFreq = 2225./freq;
t = -R:R;
% sun 2D mask initialization
sun(1:2*R+1,1:2*R+1) = 0;

% Sun 2D mask
for i = 1:2*R+1
    for j = 1:2*R+1
        if sqrt( (i-R+1)^2 + (j-R+1)^2 ) <= R
            sun(i,j) = 1;
        end
    end
end

% convolution from the Sun and vertical antenna pattern
sunShape(1:2*R+1,1:length(freq)) = 0;

for j = 1:length(freq)
    vertAntennaPattern = @(x) exp( - x.^2 / ( 2 * (DVertFreq(j)/2.355)^2 ) );
    for i = 1:2*R+1
        sunShape(i,j) = sum( sun(:,i)' .* vertAntennaPattern(t) );
    end
end

%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = 32;

for g = 0:15
    if g == 0
        maxGauss = [];
        AGauss = [];
        DGauss = [];
    else
        maxGauss(1:g) = 0;
        AGauss(1:g) = 0;
        DGauss(1:g) = 0;
    end
    
    for h = 1:g  
        %% add
        % add first gaussian
        [~, maxGauss(h)] = max(diff);
        leftBound = maxGauss(h);
        rightBound = maxGauss(h);
        
        % windowed averaging
        diff = smooth(diff,5)';

        % right bound
        while (diff(rightBound+2) - diff(rightBound+1)) < (diff(rightBound+1) - diff(rightBound))
            rightBound = rightBound + 1;
        end

        % left bound
        while (diff(leftBound-2) - diff(leftBound-1)) < (diff(leftBound-1) - diff(leftBound))
            leftBound = leftBound - 1;
        end

        xGauss = leftBound-maxGauss(h):rightBound-maxGauss(h);

        A1(1:length(xGauss)) = 1;
        A2 = xGauss.^2;
        A = [A1', A2'];
        A1 = [];
        b = log(abs(diff(leftBound:rightBound)));
        coef = ( A' * A ) \ ( A' * b' );
        AGauss(h) = exp(coef(1));
        DGauss(h) = sqrt( -1 / ( 2 * coef(2) ) );

        %% min
        % minimization with first gaussian
        %figure
        %hold on
        %axis tight

        antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
        convolution(1:length(x)) = 0;

        for j = 1:length(x)
            convolution(j) = trapz(sunShape(:,ii)' .* antennaPattern(x(j)-t));
        end
        
        gauss(1:size(data, 2)) = 0;
        for f = 1:h
            gauss = gauss + AGauss(f).*exp(-(x-maxGauss(f)+CRPIX).^2./(2.*DGauss(f).^2));
        end
        A = convolution';
        b = data(1,:,ii) - gauss;
        coef = (A'*A) \ (A'*b');

        % differences between scan and quiet Sun shape
        diff = data(1,:,ii) - coef.*convolution - gauss;

        %plot(x, data(1,:,ii))
        %plot(x, coef.*convolution + gauss, '--')
        %plot(x, diff, '-')

    end
    
    %% find
    % find new solar center
    newCRPIX = CRPIX;
    currentDiscr = discrepancy(CRPIX, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss);
    posStepDiscr = discrepancy(CRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss);
    negStepDiscr = discrepancy(CRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss);

    if posStepDiscr < currentDiscr
        while posStepDiscr < currentDiscr
            newCRPIX = newCRPIX + 1;
            currentDiscr = posStepDiscr;
            posStepDiscr = discrepancy(newCRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss);
        end
    elseif negStepDiscr < currentDiscr
        while negStepDiscr < currentDiscr
            newCRPIX = newCRPIX - 1;
            currentDiscr = negStepDiscr;
            negStepDiscr = discrepancy(newCRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss);
        end
    end

    CRPIX = newCRPIX;

    %% plot
    % plot quiet Sun shape and RATAN scan with first gaussian
    x = 1-CRPIX:1:size(data,2)-CRPIX;

    figure
    hold on
    axis tight

    antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
    convolution(1:length(x)) = 0;

    for j = 1:length(x)
        convolution(j) = trapz(sunShape(:,ii)' .* antennaPattern(x(j)-t));
    end

    gauss(1:size(data, 2)) = 0;
    for h = 1:g
        gauss = gauss + AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2));
    end
    A = convolution';
    b = data(1,:,ii) - gauss;
    coef = (A'*A) \ (A'*b');

    % differences between scan and quiet Sun shape
    diff = data(1,:,ii) - coef.*convolution;

    plot(x, data(1,:,ii))
    plot(x, coef.*convolution + gauss, '--')
    plot(x, diff - gauss, '-')
end

%% sun shape discrepancy function %%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds discrepancy from RATAN scan and quiet Sun shape
function discr = discrepancy(CRPIX, R, Dfreq, data, sunShape, ii, maxGauss, AGauss, DGauss)
    x = 1-CRPIX:1:size(data,2)-CRPIX;
    t = -R:R;
    antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
    convolution(1:length(x)) = 0;
    
    % layer with gaussians
    gauss(1:length(x)) = 0;
    for j = 1:length(maxGauss)
        gauss = gauss + AGauss(j).*exp(-(x-maxGauss(j)+CRPIX).^2./(2.*DGauss(j).^2));
    end
    
    for j = 1:length(x)
        convolution(j) = trapz(sunShape' .* antennaPattern(x(j)-t));
    end
    
    A = convolution';
    b = data(1,:,ii) - gauss;
    coef = (A'*A) \ (A'*b');
    
    diff = data(1,:,ii) - coef.*convolution - gauss;
    
    discr = sum(diff);
end