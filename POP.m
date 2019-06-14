function s=POP(im, bParams, qParams, cartImage,flagpeel)
% Inputs:
%  im       - a Velocity map image, the image is assumed to be square
%             matrix, centered and corrected for deformations, tilt, etc.
%  bParams  - a vector with the (even) beta parameter (Legendre order)
%             to be used in the fitting procedure (excluding 0th-order
%             which is always included), e.g. [2 4] for beta_2 & beta_4.
%  qParams  - a vector that states the quadrants of the image to include
%             in the analysis, e.g. [2 3] will treat only the left side
%             of the image. This is sometimes needed if there is a burned
%             spot, or for images that are only symmetric in half of the
%             plane, as in two color (w,2w) experiments
%  cartImage- the kind of Cartesian image to return, acceptable value are
%             'sim' for the simulated (fit) image
%             'exp' for the experimental image
%              0   for no Cartesian image (speeds up run time!)
% flagpeel -   1 for onion peeling in R + Legendre decomposition in theta,
%              0 for only Legendre decomposition
%
% Outputs:
%  s            - A Matlab structure that contains the following:
%  s.iraraw     - a 2d triangular array (polar) of the raw data folded
%                 to (0,pi/2) - G(R,alpha)
%  s.iradecon   - a 2d triangular array (polar) of the simulated
%                 deconvolved data folded to (0,pi/2) - g_fit(r;R,alpha)
%  s.iradeconExp- a 2d triangular array (polar) of the experimental
%                 deconvolved data folded to (0,pi/2) - g_fit(r;R,alpha)
%  s.sqp        - a 2d rectangular array (polar) derived from iraraw
%  s.PESR       - a 1d vector of the radial projection of the raw data
%  s.PESId      - a 1d vector of the radial projection of the (simulated)
%                 deconvolved data
%  s.PESIdExp   - a 1d vector of the radial projection of
%                 the (experimental) deconvolved data
%  s.Betas      - a matrix containing each beta parameters (as specified in
%                 bParams, e.g. if bParams=[2 4] then  s.Betas(1,:) is beta_0,
%                 s.Betas(2,:) gives b2 and s.Betas(3,:) gives b4
%
% And depending on the option selected for cartImage:
%  s.simImage    - a 2d reconstruction in Cartesian
%                  coordinates of the simulated image
%  s.expImage    - a 2d reconstruction in Cartesian
%                  coordinates of the experimental image
%  s.sqpdecon    - a 2d rectangular array (polar) derived from iradecon
%  s.sqpdeconExp - a 2d rectangular array (polar) derived from iradeconExp
%
%   Example: see the script pop_example.m 

%   Ver 1.7 (2018-10-01)
%   Adi Natan (natan@stanford.edu)
%
%% defaults
if (nargin < 5);                                               flagpeel = 1; end
if (nargin < 4);                                cartImage = 0; flagpeel = 1; end
if (nargin < 3);                  qParams=1:4 ; cartImage = 0; flagpeel = 1; end
if (nargin < 2);  bParams=[2 4];  qParams=1:4 ; cartImage = 0; flagpeel = 1; end
if (nargin < 1);  load testimg ; bParams=[2 4];  qParams=1:4 ; cartImage = 'sim'; flagpeel = 1; im = im'; end


% Check that the beta Parameters are in the scope of the code
if any(mod(bParams,2)) || any(bParams<=0) || any(bParams>42)
    error('Only even positive beta parameters of <=42 orders supported! Beyond that order there are floating point accuracy errors and vpa treatment is needed');
end

% Check that the images size in the scope of the code
if max(size(im))>4156
    error('incompatible image size, the code supports images up to 4156x4156 pixels');
end

if min(size(im))<15
    error('image size too small to obtain meaningful analysis');
end


%% Here we go:
x0=ceil(size(im,1)/2); y0=ceil(size(im,2)/2); % Assuming image is centered
L = min([x0,y0]);%,(size(im) - [x0,y0])]);%set radius of square matrix around center
RR=(0:size(im)/2);
PPR=(floor(0.5*pi*(RR+1))-1); % # of pixels per radius
AngleInc = (0.5*pi./PPR'); % angle increment per radius
AngleInc(1)=0; % avoid inf at origin

iraraw=tripolar(im,qParams); % go to polar coordinates via tripolar
PESR=nansum(iraraw,2)./nansum(iraraw>0,2).*(1:size(iraraw,1));


if flagpeel
    % load radial basis set  (radial projection of spherical delta functions)
    matObjLut = matfile('lut.mat');
    lut=matObjLut.lut(1:L,1:L); % this loads only the relevant size lut
    
    % Do the Legendre decomposition
    [betas, iradecon, iradeconExp, PESId, PESIdExp]=LDP(iraraw, bParams,lut);
else % Legendre decomposition without peeling
    [betas, iradecon, iradeconExp, PESId, PESIdExp]=LDP(iraraw, bParams);
end


PESId(~isfinite(PESId))=0;
iradecon(~isfinite(iradecon))=0;
iradeconExp(~isfinite(iradeconExp))=0;

%  Create a square polar projection from the triangular one
sqp=t2s(iraraw);
sqpd=t2s(iradecon);
sqpde=t2s(iradeconExp);

% 2d transform to Cartesian coordinates of the simulated\exp image
if strcmp(cartImage, 'sim')
    zz=t2f(iradecon,L,PPR,AngleInc,im);
    s=struct('iraraw',iraraw,'iradecon',iradecon, 'iradeconExp', iradeconExp,...
        'sqp',sqp,'sqpd',sqpd,'sqpde',sqpde,'PESR',PESR,'PESId',PESId,'PESIdExp',PESIdExp,...
        'Betas',betas,'simImage',zz);
    
elseif (strcmp(cartImage,'exp'))
    zz=t2f(iradeconExp,L,PPR,AngleInc,im);
    s=struct('iraraw',iraraw,'iradecon',iradecon, 'iradeconExp', iradeconExp,...
        'sqp',sqp,'sqpd',sqpd,'sqpde',sqpde,'PESR',PESR,'PESId',PESId,'PESIdExp',PESIdExp,'Betas',betas,'expImage',zz);
    
elseif cartImage==0
    s=struct('iraraw',iraraw,'iradecon',iradecon,'iradeconExp', iradeconExp,...
        'sqp',sqp,'sqpd',sqpd,'sqpde',sqpde,'PESR',PESR,'PESId',PESId,'PESIdExp',PESIdExp, 'Betas',betas);
end

% demo mode - in case there's no input to the function, use test img
if (nargin < 1)
    figure;
    imagesc([im(:,1:end/2)./max(max(im(:,1:end/2))) s.simImage(:,end/2+1:end)./max(max( s.simImage(:,end/2+1:end)))]); axis square
    title('raw vs onion peeled');
    
    figure('Position',[0 0 1000 600]);
    r=1:numel(s.PESId);
    subplot(2,3,1);imagesc(im); title('raw image');xlabel('x-pixels');ylabel('y-pixels');
    subplot(2,3,2);imagesc(s.iraraw); title('iraraw');xlabel('folded angle (PPR)');ylabel('radius');
    subplot(2,3,3);imagesc([0 pi/2],[0 size(s.sqp,1)],s.sqp); title('sqp');xlabel('angle [rad]');ylabel('radius');
    
    subplot(2,3,4);imagesc(s.simImage); title('sim Image');xlabel('x-pixels');ylabel('y-pixels');
    subplot(2,3,5);imagesc(s.iradecon);title('iradecon');xlabel('folded angle (PPR)');ylabel('radius');
    subplot(2,3,6);imagesc([0 pi/2],[0 size(s.sqpd,1)],s.sqpd); title('sqpd');xlabel('angle');ylabel('radius');
    
    sb=size(s.Betas,1);
    figure('Position',[0 0 250*(sb+1) 250]);
    subplot(1,sb+1,1); plot(r,s.PESId);title('PESId');xlabel('radius');ylabel('intensity');
    for nsb=1:sb
        subplot(1,sb+1,1+nsb);  plot(r,s.Betas(nsb,:).*(abs(s.PESId)>0.1*max(abs(s.PESId)))); title(['\beta_{' num2str(nsb*2-2) '}']);xlabel('radius');ylabel('intensity');
    end
end

function ira=tripolar(im, qParams)
% This function accepts a centered image (im), and quadrant parameters (1
% to 4) and produces a polar representation that transforms each Cartesian
% quadrant to the correct polar pixel size.
% The code starts with a Cartesian image of square pixel size dx*dy=1, and
% outputs a polar image of polar pixels that has a similar amount of
% information in polar coordinates dr*dtheta~=1.  The signal in each polar
% pixel is determined via its fractional overlap with the four surrounding
% Cartesian pixels. The result is a triangular polar representation,
% because the number of polar pixels per radius per angle increase with
% radius, i.e. the larger the radius the higher the angular resolution.
% The quadrants are then averaged because of the cylindrical symmetry in VMI.
% The code supports NaN valued pixels for masking.

x0=ceil(size(im,1)/2); y0=ceil(size(im,2)/2); % Assuming image is centered
RR=(0:size(im)/2);
PPR=(floor(0.5*pi*(RR+1))-1); % calc the  # of pixels per radius
AngleInc = (0.5*pi./PPR'); % angle increment per radius
AngleInc(1)=0; % avoid inf at origin

%set radius of square matrix around center
L = min([x0,y0,(size(im) - [x0,y0])]);

%create  image quadrants
if mod(size(im,1),2)==1 %case for odd pixel image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0:y0+L-1);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0:x0+L-1   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0:x0+L-1   ,y0:y0+L-1);
else %case for even pixel image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0+1:y0+L);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0+1:x0+L   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0+1:x0+L   ,y0+1:y0+L);
end

%add image quadrants.  NaN's are not counted.  If more than one
%quadrant contributes, values are averaged. This assumes that each quadrant
%has the same response or sensitivity. 
% NOTE: zero values are assumed to be valid values, this is to allow
% analysis of difference images. For masking use NaNs. 
a4=nansum(Q(:,:,qParams),3)./sum( ~isnan(Q(:,:,qParams)).*(Q(:,:,qParams)>-inf),3);

ira = NaN(L-2,PPR(L)); % initialize the  matrix
ira(1,1) = a4(1,1);    % origin pixel remains the same

% creating the 2d triangular array polar image
for r=2:L-2
    npr=PPR(r); % determine # polar pix in radius
    angincr=AngleInc(r);    % the angular increment per radius
    qp=0:npr;
    xp=r*sin(angincr*qp)+1;  % polar x-coordinate
    yp=r*cos(angincr*qp)+1;  % polar y-coordinate
    
    % define scale fractional weight of cart pixels in polar pixels
    xc=round(xp);yc=round(yp);
    xd=1-abs(xc-xp);
    yd=1-abs(yc-yp);
    
    a4r=[ (a4(xc+(yc-1)*L));
        0 (a4(xc(2:npr)+(yc(2:npr)-1+(-1).^(yp(2:npr)<yc(2:npr)))*L)) 0;
        0 (a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+yc(2:npr)*L)) 0;
        0 (a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+L*(yc(2:npr)+(-1).^(yp(2:npr)<yc(2:npr))))) 0];
    
    c=[xd.*yd;
        0 xd(2:npr).*(1-yd(2:npr)) 0;
        0 (1-xd(2:npr)).*yd(2:npr) 0;
        0 (1-xd(2:npr)).*(1-yd(2:npr)) 0];
    
    c=bsxfun(@rdivide,c.*(~isnan(a4r)),sum(c.*(~isnan(a4r))));
    
    ira(r,1:npr+1) = nansum(c.*a4r);
    % only when all four surrounding Cartesian pixels are NaN assign NaN
    ira(r,all(isnan(c.*a4r)))=NaN;
    
end

function [betas, iradecon, iradeconExp, PESId, PESIdExp]=LDP(ira, bParams,lut)
% This function accepts the triangular polar representation (ira) from the
% previous function tripolar, as well as beta parameters (bParams) 
% for Legendre decomposition in angle and a basis set (lut) for onion peeling
% The outputs are detailed in the header of this code

% defaults
if (nargin < 2)
    bParams=[2 4];
    flagpeel=0;
end

if exist('lut','var')
    flagpeel=1;
else
    flagpeel=0;
end

PPR=(floor(0.5*pi*(1:size(ira,1) ))-1); % calc the  # of pixels per radius
AngleInc = 0.5*pi./PPR'; % angle increment per radius
AngleInc(1)=0;
L=size(ira,1);
betas = zeros(numel(bParams)+1,L);
PESId = zeros(1,L);

for r =  L-2:-1:2
    
    npr=PPR(r); % # of polar pixels in radius -1
    qp=0:npr;
    
    % for very small radii we reduce the maximal beta value
    % allowed to avoid divergences at the origin
    if npr/2 <= max(bParams)
        bParams=2:2:npr/2;
    end
    
    y = ira(r,qp+1); % assign row of all pixels in radius rp
    % this may have NaN values, one way to treat such a case is to 
    % interpolate \ extrapolate NaN values using splines or AR methods. 
    % We will ignore this solution for now and keep the NaNs as is.
    
    B = zeros(1,numel(bParams)+1);
    % this is beta0, including if there are NaN values
    B(1)=nansum(y)/sum(~isnan(y))*(npr+1);
    
    % one fit coefficient for each B param
    fitCoefs = ones(numel(bParams) + 1, sum(~isnan(y)));
    for ii=(1:numel(bParams))
        % assign relevant Legendre polynomials to fitCoef
        fitCoefs(ii+1,:) = leg(bParams(ii), cos(AngleInc(r)*qp(~isnan(y))));
        B(ii+1) = y(~isnan(y))*fitCoefs(ii+1,:)'; % fill B vector
    end
    
    A = fitCoefs * fitCoefs'; % matrix A for least square fitting
    A(1,1)=npr+1;
    lastwarn(''); % reset last warning
    
    Ain=A\eye(size(A)); % instead of inv(A)
    
    [~, msgid] = lastwarn; % capture warning in case Ain is close to singular
    if strcmp(msgid,'MATLAB:nearlySingularMatrix')
        % switch to Moore-Penrose pseudoinverse in case of nearly singular
        Ain =  pinv(A);
    end
    
    Beta = zeros(1,numel(bParams)+1);
    Beta(1)=B*Ain(:,1);
    
    betas(1,r)=(Beta(1));   % Beta0 is just the Intensity Factor
    
    
    for ii=1:numel(bParams)
        Beta(ii+1)=B*Ain(:,ii+1)/Beta(1); % generate beta matrix
        %  Beta(ii+1)=B*Ain(:,ii+1); % generate beta matrix
        betas(ii+1,r) = Beta(ii+1); % copy for output betas
    end
    
    %%%%%%%
    
    if 0==Beta(1)
        PESId(r)=0;
        continue;
    end
    % generate matrices for alpha; R/rp * cos(alpha); and the basis set
    % scaled by pixels per radius
    alphaMat = (AngleInc(1:r) * (0:PPR(r)));
    rrpCosAlphaMat = repmat((1:r)'/r, 1, PPR(r)+1).*cos(alphaMat);
    if flagpeel
        itMat = repmat((lut(r+1,2:r+1)*(npr+1)./ (PPR(1:r)+1))', 1, PPR(r)+1);
    end
    bContrib = ones(r, PPR(r)+1);
    
    for ii=1:numel(bParams) %  add each beta contribution
        bContrib = bContrib + Beta(ii+1)*leg(bParams(ii), rrpCosAlphaMat);
    end
    
    % generate the simulated image for this radius
    if flagpeel
        factMat = Beta(1).*itMat.*bContrib;
    else
        factMat = Beta(1).*bContrib;
    end
    
    % save the simulated data
    PESId(r) = PESId(r) + sqrt(r)*factMat(end,:)*sin(alphaMat(end,:))';
    iradecon(r,1:PPR(r)+1)=factMat(end,1:PPR(r)+1);%/sqrt(r);
    % save the experimental data
    PESIdExp(r) = PESId(r) + sqrt(r)*ira(r,1:PPR(r)+1)*sin(alphaMat(end,:))';
    iradeconExp(r,1:PPR(r)+1)=ira(r,1:PPR(r)+1);%/sqrt(r);
    %% subtract away the simulated image
    if flagpeel
        % ira(1:r,1:PPR(r)+1) = ira(1:r,1:PPR(r)+1)-factMat;
%        ira(1:r,qp+1) = max(ira(1:r,qp+1)-factMat,zeros(r,numel(qp))); % subtract fact from ira unless smaller than zero
        ira(1:r,qp+1) = ira(1:r,qp+1)-factMat; % subtract fact from ira unless smaller than zero
    end
end


function squmat=t2s(trimat)
% transfer from triangular to square polar matrix 
% this is just done for visualizing the polar data
PPR=(floor(0.5*pi*((0:size(trimat,1))+1))-1); % calc the  # of pixels per radius in quadrant
[Ms,Ns] = size(trimat);
% create a matrix that "squares" the triangular polar matrix
squmat=zeros(Ms,Ns);
squmat(1,:)=ones(1,Ns).*trimat(1,1);

for k=2:Ms
    try % in case data is all zero, then interp1 wont work
        %first map pixels in given range to new range
        smap=interp1(1:PPR(k) , single(~isnan(trimat(k,1:PPR(k)))) , linspace(1,PPR(k),Ns) ,'nearest');
        % interp only on non NaN pixels
        nr=~isnan(trimat(k,1:PPR(k)));
        squmat(k,:)=interp1(find(nr),trimat(k,nr),linspace(1,PPR(k),Ns),'nearest');
        squmat(k,~smap)=NaN;
        
    catch
    end
end

function fullcart=t2f(trimat,L,PPR,AngleInc,im)
% transfer from triangular polar matrix to full Cartesian matrix
x = []; y = []; z = [];
for r = 1:L-5
    qp = 0:PPR(r);
    x = [x r*cos((qp)*AngleInc(r))];
    y = [y r*sin((qp)*AngleInc(r))];
    z = [z trimat(r,qp+1)];
end

try % Making sure the interpolation function is Matlab version independent
    F = scatteredInterpolant(double(x(:)),double(y(:)),double(z(:)),'natural');
catch
    disp('for older matlab versions use F = TriScatteredInterp(...)');
end

[xx,yy] = meshgrid(1:L,1:L);
zz = F(xx,yy);
%zz = max(zz,zeros(size(zz)));

m4=@(m)[rot90(m,2), flipud(m); fliplr(m),m]; % fold back from quad to full
fullcart=m4(zz);

if all(mod(size(im),2)) %for the case of odd # of pixels in image
    fullcart(end/2,:)=[];
    fullcart(:,end/2)=[];
end

% masking the central pixels where noise accumulates, and anything beyond L
[xm,ym] = meshgrid(1:size(fullcart,1));
innermask=1-exp(-(xm-size(fullcart,1)/2).^2/10-(ym-size(fullcart,2)/2).^2/10);
outermask=double((xm-size(fullcart,1)/2).^2+(ym-size(fullcart,1)/2).^2 <L.^2);
outermask(outermask==0)=NaN;
fullcart=fullcart.*outermask.*innermask;

function p=leg(m,x)
%  This function returns Legendre polynomial P_m(x) where m is the degree
%  of polynomial and X is the variable.
%  x2 is used for improved performance by minimizing the # of operations.
switch m
    case 0
        p=ones(size(x));
        return
    case 1
        p=x;
        return
    case 2
        p=(3*x.*x -1)/2;
        return
    case 4
        x2=x.*x;
        p = ((35.*x2-30).*x2+3)/8;
        return
    case 6
        x2=x.*x;
        p = (((231.*x2-315).*x2+105).*x2-5)/16;
        return
    case 8
        x2=x.*x;
        p = ((((6435.*x2-12012).*x2+6930).*x2-1260).*x2+35)/128;
        return
    case 10
        x2=x.*x;
        p = (((((46189.*x2-109395).*x2+90090).*x2-30030).*x2+3465).*x2-63)/256;
        return
    case 12
        x2=x.*x;
        p = ((((((676039.*x2-1939938).*x2+2078505).*x2-1021020).*x2+225225).*x2-18018).*x2+231)/1024;
        return
    case 14
        x2=x.*x;
        p = (((((((5014575.*x2-16900975).*x2+22309287).*x2-14549535).*x2+4849845).*x2-765765).*x2+45045).*x2-429)/2048;
        return
    case 16
        x2=x.*x;
        p = ((((((((300540195.*x2-1163381400).*x2+1825305300).*x2-1487285800).*x2+669278610).*x2-162954792).*x2+19399380).*x2-875160).*x2+6435)/32768;
        return
    case 18
        x2=x.*x;
        p = (((((((((2268783825.*x2-9917826435).*x2+18032411700).*x2-17644617900).*x2+10039179150).*x2-3346393050).*x2+624660036).*x2-58198140).*x2+2078505).*x2-12155)/65536;
        return
    case 20
        x2=x.*x;
        p = ((((((((((34461632205.*x2-167890003050).*x2+347123925225).*x2-396713057400).*x2+273491577450).*x2-116454478140).*x2+30117537450).*x2-4461857400).*x2+334639305).*x2-9699690).*x2+46189)/262144;
        return
    case 22
        x2=x.*x;
        p = (((((((((((263012370465.*x2-1412926920405).*x2+3273855059475).*x2-4281195077775).*x2+3471239252250).*x2-1805044411170).*x2+601681470390).*x2-124772655150).*x2+15058768725).*x2-929553625).*x2+22309287).*x2-88179)/524288;
        return
    case 24
        x2=x.*x;
        p = ((((((((((((8061900920775.*x2-47342226683700).*x2+121511715154830).*x2-178970743251300).*x2+166966608033225).*x2-102748681866600).*x2+42117702927300).*x2-11345993441640).*x2+1933976154825).*x2-194090796900).*x2+10039179150).*x2-202811700).*x2+676039)/4194304;
        return
    case 26
        x2=x.*x;
        p= (((((((((((((61989816618513 .*x2-395033145117975) .*x2+1112542327066950) .*x2-1822675727322450) .*x2+1923935489951475) .*x2-1369126185872445) .*x2+667866432132900) .*x2-222622144044300) .*x2+49638721307175) .*x2-7091245901025) .*x2+601681470390) .*x2-26466926850) .*x2+456326325) .*x2-1300075)/8388608;
        return
    case 28
        x2=x.*x;
        p= ((((((((((((((956086325095055 .*x2-6570920561562378) .*x2+20146690401016725) .*x2-36343049350853700) .*x2+42832879592077575) .*x2-34630838819126550) .*x2+19624141997505045) .*x2-7823578204985400) .*x2+2170565904431925) .*x2-408140597414550) .*x2+49638721307175) .*x2-3610088822340) .*x2+136745788725) .*x2-2035917450) .*x2+5014575)/33554432;
        return
    case 30
        x2=x.*x;
        p= (((((((((((((((7391536347803839 .*x2)-54496920530418135 .*x2)+180700315442965395 .*x2)-355924863751295475 .*x2)+463373879223384675 .*x2)-419762220002360235 .*x2)+271274904083157975 .*x2)-126155198555389575 .*x2)+42051732851796525 .*x2)-9888133564634325 .*x2)+1591748329916745 .*x2)-166966608033225 .*x2)+10529425731825 .*x2)-347123925225 .*x2)+4508102925 .*x2-9694845)/67108864;
        return
    case 32
        x2=x.*x;
        p= ((((((((((((((((916312070471295267 .*x2)-7214139475456546864 .*x2)+25722546490357359720 .*x2)-54932895894661480080 .*x2)+78303470025285004500 .*x2)-78588209916286040880 .*x2)+57087661920320991960 .*x2)-30382789257313693200 .*x2)+11858588664206620050 .*x2)-3364138628143722000 .*x2)+680303589246841560 .*x2)-94926082220489520 .*x2)+8682263617727700 .*x2)-479493848710800 .*x2)+13884957009000 .*x2)-158685222960 .*x2+300540195)/2147483648;
        return
    case 34
        x2=x.*x;
        p= (((((((((((((((((7113260368810144185 .*x2)-59560284580634192355 .*x2)+227245393476881226216 .*x2)-523025111970599647640 .*x2)+810260214446256831180 .*x2)-892659558288249051300 .*x2)+720391924232622041400 .*x2)-432235154539573224840 .*x2)+193690281515374794150 .*x2)-64563427171791598050 .*x2)+15811451552275493400 .*x2)-2783060137827988200 .*x2)+340151794623420780 .*x2)-27382523717448900 .*x2)+1335732864265800 .*x2)-34249560622200 .*x2)+347123925225 .*x2-583401555)/4294967296;
        return
    case 36
        x2=x.*x;
        p= ((((((((((((((((((110628135069209194801 .*x2)-981629930895799897530 .*x2)+3990539066902490887785 .*x2)-9847300383998186469360 .*x2)+16475291027073888900660 .*x2)-19770349232488666680792 .*x2)+17555637979668898008900 .*x2)-11732097051788416102800 .*x2)+5943233374919131841550 .*x2)-2281241093403303131100 .*x2)+658546957152274300110 .*x2)-140865659283908941200 .*x2)+21800637746319240900 .*x2)-2354897039700605400 .*x2)+168206931407186100 .*x2)-7302006324653040 .*x2)+166966608033225 .*x2)-1511010027450 .*x2+2268783825)/17179869184;
        return
    case 38
        x2=x.*x;
        p= (((((((((((((((((((861577581086657669325 .*x2)-8075853860052271220473 .*x2)+34847862546800896362315 .*x2)-91782398538757290419055 .*x2)+164942281431969623361780 .*x2)-214178783351960555708580 .*x2)+207588666941131000148316 .*x2)-152984845251400396934700 .*x2)+86524215756939568758150 .*x2)-37640478041154501663150 .*x2)+12546826013718167221050 .*x2)-3172998975370048900530 .*x2)+598679051956613000100 .*x2)-82171634582280215700 .*x2)+7905725776137746700 .*x2)-504620794221558300 .*x2)+19624141997505045 .*x2)-402684172315425 .*x2)+3273855059475 .*x2-4418157975)/34359738368;
        return
    case 40
        x2=x.*x;
        p=((((((((((((((((((((26876802183334044115405 .*x2)-265365894974690562152100 .*x2)+1211378079007840683070950 .*x2)-3391858621221953912598660 .*x2)+6516550296251767619752905 .*x2)-9104813935044723209570256 .*x2)+9566652323054238154983240 .*x2)-7710436200670580005508880 .*x2)+4819022625419112503443050 .*x2)-2345767627188139419665400 .*x2)+888315281771246239250340 .*x2)-260061484647976556945400 .*x2)+58171647881784229843050 .*x2)-9763073770369381232400 .*x2)+1197358103913226000200 .*x2)-103301483474866556880 .*x2)+5929294332103310025 .*x2)-207785032914759300 .*x2)+3847870979902950 .*x2)-28258538408100 .*x2+34461632205)/274877906944;
        return
    case 42
        x2=x.*x;
        p=(((((((((((((((((((((209863810776486386280915 .*x2)-2177020976850057573347805 .*x2)+10481952851500277205007950 .*x2)-31092037361201244198821050 .*x2)+63597349147911635861224875 .*x2)-95141634325275807248392413 .*x2)+107740298231362557979914696 .*x2)-94299858612963204670549080 .*x2)+64574903180616107546136870 .*x2)-34804052294693590302644250 .*x2)+14778336051285278343892020 .*x2)-4926112017095092781297340 .*x2)+1278635632852551404981550 .*x2)-255060302250900084696450 .*x2)+38354932669308283413000 .*x2)-4230665300493398534040 .*x2)+329273478576137150055 .*x2)-17090318957238952425 .*x2)+542549808166315950 .*x2)-9113378636612250 .*x2)+60755857577415 .*x2-67282234305)/549755813888;
        return
end

function y = nansum(x,dim)
% Sum, ignoring NaNs.
x(isnan(x)) = 0;
if nargin == 1 % let sum figure out which dimension to work along
    y = sum(x);
else           % work along the explicitly given dimension
    y = sum(x,dim);
end
