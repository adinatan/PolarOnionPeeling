function [trimat, varargout]=beta2cart(betas,bParams)
% this function transform from the beta coeficcients that are associated 
% with Legendre decomposition to a tripolar and cartesian representation
% in 2D and 3D.
%
% Inputs:
% betas - a matrix of beta parameters betas(bParams,r), 
% bParams=0,2,4... the Legendre order associated with each beta
%
%Outputs:
% trimat - is a tripolar representation betweeb 0 and pi/2
% varagout(1) - a 2D cartesian representation
% varagout(2) - a 3D cartesian representation
%
% for proper conversion one needs at least beta0, becuase by defenition
% the Legendre decomposition is I(r)=beta0(r)*(1+sigma(beta(r)_2n*leg(2n)))
%
% % example:
% % generate data
% r=linspace(0,10,64); % radius vec
% b0=exp(-(r-5).^2/5); % beta 0
% b2=2*b0.*besselj(2,r); % beta 2
% b4=b0.*besselj(4,r); % beta 4
% beta_in_vec=[b0(:)';b2(:)';b4(:)']; % create betas input vector
% 
% [s1, s2, s3]=beta2cart(beta_in_vec); % make it happen
% xx=[-fliplr(r) r(2:end)]; % from radial to carteisan range
% 
% % show all outputs
% subplot(2,2,1)
% plot(r,beta_in_vec'); axis square;
% xlabel('radius'); ylabel('\beta amplitude');
% title('beta input vector');
% 
% subplot(2,2,2);
% imagesc(linspace(0,pi/2,size(s1,2)),r,s1); axis square;
% xlabel('Angle [rad]') ; ylabel('radius');
% title('s(1)');
% 
% subplot(2,2,3);
% imagesc(xx,xx,s2'); axis square;
% xlabel('x');ylabel('y')
% title('s(2)');
% 
% subplot(2,2,4);
% s3(s3<0.5*max(s3(:)))=nan; % take out small values for visualization
% [X, Y ,Z]=meshgrid(xx);
% h = slice(X,Y,Z,s3, [], [], xx);
% set(h, 'EdgeColor','none', 'FaceColor','interp');
% % set view
% axis tight; daspect([1 1 1])
% view([145,-15])
% 
% 
% xlim([xx(round(numel(xx)*0.15)) xx(round(numel(xx)*0.85)) ]);
% ylim([xx(round(numel(xx)*0.15)) xx(round(numel(xx)*0.85)) ]);
% zlim([xx(round(numel(xx)*0.15)) xx(round(numel(xx)*0.85)) ]);
% 
% % normalize s3 for alpha (transparency) map
% ns3=(s3(:)-min(s3(:)))./(max(s3(:))-min(s3(:)));
% ns3=reshape(ns3,size(s3));
% % set alpha map
% for n=1:numel(h)
%     set(h(n),'AlphaData',ns3(:,:,n),'FaceAlpha','interp')
% end
% title('s(3)');
% __________________________________
%   Adi Natan (natan@stanford.edu)
%   Ver 1.2 , Date: Jul 11th 2019

if (nargin < 2)
    bParams=2:2:size(betas,1)*2-2; % ignoring beta0, consitant with LD code
end

L = size(betas,2);
PPR=(floor(0.5*pi*((0:L)+1))-1); % # of pixels per radius
AngleInc = single(0.5*pi./PPR'); % angle increment per radius
AngleInc(1)=0; % avoid inf at origin

trimat=NaN(L,PPR(L));

for r = 2:L
    alphaMat = (AngleInc(1:r) * (0:PPR(r)));
    rrpCosAlphaMat = repmat((1:r)'/r, 1, PPR(r)+1).*cos(alphaMat);
    clear beta_contr
    
    for ii=1:numel(bParams) %  add each beta contribution
        beta_contr(:,:,ii)= betas(ii+1,r)*leg(bParams(ii), rrpCosAlphaMat);
    end
    
    factMat = betas(1,r).*(ones(r, PPR(r)+1) + sum(beta_contr,3));
    % note that the /sqrt(r) is omitted because the function is suppose
    % to work such as LD2(beta2cart(betas))=betas;
    trimat(r,1:PPR(r)+1)=factMat(end,1:PPR(r)+1);% /sqrt(r);
end

if max(nargout,1)>1
    nout = max(nargout,1) - 1;
    if nout
        %% 2D- cartesian:
        x = []; y = []; z = [];
        for r = 1:L
            qp = 0:PPR(r);
            x = [x r*cos(qp*AngleInc(r))];
            y = [y r*sin(qp*AngleInc(r))];
            z = [z trimat(r,qp+1)];
        end
        
        try % Making sure the interpolation function is Matlab version independent
            F = scatteredInterpolant(double(x(:)),double(y(:)),double(z(:)),'natural');
        catch
            F = TriScatteredInterp(double(x(:)),double(y(:)),double(z(:)),'natural');
        end
        
        [xx,yy] = meshgrid(0:L-1,0:L-1);
        zz = F(xx,yy);
        zz = max(zz,zeros(size(zz)));
        zz((xx.^2+yy.^2)>(L-1).^2)=0; % nothing outside L
        sb=size(betas,2);
        fullcart = zeros(2*sb-1);
        fullcart(1:(size(betas,2)),1:(size(betas,2)))= zz(end:-1:1,end:-1:1);
        fullcart((size(betas,2)):end,1:(size(betas,2)))= zz(:,end:-1:1);
        fullcart(1:(size(betas,2)),(size(betas,2)):end)= zz(end:-1:1,:);
        fullcart((size(betas,2)):end,(size(betas,2)):end)= zz;
        varargout{1}=fullcart;%.*mask;
        
        %% from 2D cartesian to 3D
        I=zz'; % make cos^2 now in vertical direction in 3d, note that originally for Legendre decompostion it was in the horizontal direction
        
        s=size(I);
        % V2 is the volume 3D matrix
        V2 = zeros(2*s(1)-1,2*s(2)-1,s(1));
        
        [i1,i2,i3] = ndgrid(1:2*s(1)-1,1:2*s(2)-1,1:s(1));
        % it will be assumed that input 1:N span elements from 0 to N-1
        % output V spans -(N-1):N-1 along x and y
        x = i1-s(1); % -(N-1):N-1
        y = i2-s(2); % -(N-1):N-1
        z = i3-1; %      0:N-1
        % input array I is in xz plane, rotated along z axis, geometrically speaking
        
        % identify the cylindrical coordinates of each voxel [i1,i2,i3]
        [~,r_out,z_out] = cart2pol(x,y,z); % theta is redundant; z_out===z
        
        % identify the coordinates of each input pixel with the above
        [j1,j2] = meshgrid(1:s(1),1:s(2));
        r_in = j1-1; % Cartesian input x <-> cylindrical output r
        z_in = j2-1; % Cartesian input y <-> cylindrical output z
        % note that j1 and j2 are swapped with respect to x and y
        % but this is what interp2 will expect later
        
        % interpolate each voxel based on r and z
        %method = 'nearest';  
        %method = 'spline';  
        method = 'cubic';  
        
        V2(:) = interp2(r_in,z_in,I,...
            r_out(:),z_out(:),method,...
            0); % extrapolation value, otherwise NaNs appear outside
        V2 = cat(3, V2(:,:,end:-1:2),V2);
        varargout{2} = permute(V2,[3 2 1]); % back to cos^2 in the horiz plane
    end
end


function [betas varargout]=LDSD(vec, bParams)
warning('off','all')

% Apply Legendre decomposition for a single data vector spanning [0,2*pi] 
% Inputs:
%  vec      - a vector that contain information that spans 0 to 2*pi
%  bParams  - a vector with the (even) beta parameter (legenadre order)
%             to be used in the fitting procedure (excluding 0th-order
%             which is always included), e.g. [2 4] for beta_2 & beta_4.
%
% Outputs:
%  betas      - a vector containing each beta_n such that positive
%               orders are normalized by beta_0 (intensity of vec),
%               for example: betas(1) is beta_0, betas(2) is beta_2/beta_0,
%               betas(3) is beta4/beta_0, ...
%
% reco_vec    - if a second output is asked,  the reconstructed vector 
%               using the betas is given 
%   Ver 1 (2019-11-14)
%   Adi Natan (natan@stanford.edu)
%
%
% % example:
% 
% % generate random vec up to beta 4 + noise
% N=64; % # of bins in angle from 0 to 2*pi
% vec=1000*(rand(1)-0.5)+...
%     (rand(1)-0.5)*cos(linspace(0,2*pi,N)).^2+...
%     (rand(1)-0.5)*cos(linspace(0,2*pi,N)).^4+...
%     0.1*(rand(1)-0.5)*rand(1,N);
% 
% %apply Legendre decomposition up to a higher order and see that indeed
% %higher order dont contribute
% [betas, reco_vec]=LDSD(vec, 2:2:6);
% 
% plot(linspace(0,2*pi,N),vec,'x');hold on; plot(linspace(0,2*pi,N),reco_vec); legend('vec','reco vec');
% for n=1:numel(betas)
%     if n>1
%         S{n}=['\beta_' num2str(2*n-2) ' = ' num2str(betas(n)*betas(1))]
%     else
%         S{n}=['\beta_' num2str(2*n-2) ' = ' num2str( betas(n))]
%     end
% end
% text(1.9*pi,betas(1),S)


%% defaults

if (nargin < 2);  bParams=[2 4]; end

% Check that the beta Parameters are in the scope of the code
if any(mod(bParams,2)) || any(bParams<=0) || any(bParams>42)
    error('Only even positive beta parameters of <=42 orders supported! Beyond that order there are floating point accuracy errors and vpa treatment is needed');
end

PPR=numel(vec)-1;
AngleInc = 2*pi./PPR'; % angle increment per radius
betas = zeros(numel(bParams)+1,1);


npr=PPR;
qp=0:npr;

% for very small radii we reduce the maximal beta value
% allowed to avoid divergences at the origin
if npr/2 <= max(bParams)
    bParams=2:2:npr/2;
end

y = vec;%ira(r,qp+1); % assign row of all pixels in radius rp
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
    fitCoefs(ii+1,:) = leg(bParams(ii), cos(AngleInc*qp(~isnan(y))));
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

betas(1)=(Beta(1));   % Beta0 is just the Intensity Factor


for ii=1:numel(bParams)
    Beta(ii+1)=B*Ain(:,ii+1)/Beta(1); % generate beta matrix
    %  Beta(ii+1)=B*Ain(:,ii+1); % generate beta matrix
    betas(ii+1) = Beta(ii+1); % copy for output betas
end



%% reconstruct signal
if nargout>1
alphaMat = (AngleInc * (0:PPR));
        rrpCosAlphaMat = ones(1,PPR+1).*cos(alphaMat);
        clear beta_contr
        
        for ii=1:numel(bParams) %  add each beta contribution
            beta_contr(:,ii)= betas(ii+1)*leg(bParams(ii), rrpCosAlphaMat);
        end
        
        factMat = betas(1).*(ones(1, PPR+1) + sum(beta_contr,2));
        % note that the /sqrt(r) is omitted because the function is suppose
        % to work such as LD2(beta2cart(betas))=betas;
        reco_vec(1:PPR+1)=factMat(1:PPR+1,1);% /sqrt(r);
        varargout{1}=reco_vec;
end
    

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
if any(mod(bParams,2)) || any(bParams<=0)  
    error('Only even positive beta parameters supported!');
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
  
        P = legendre(m,x);
        p = squeeze(P(1,:,:));  


function y = nansum(x,dim)
% Sum, ignoring NaNs.
x(isnan(x)) = 0;
if nargin == 1 % let sum figure out which dimension to work along
    y = sum(x);
else           % work along the explicitly given dimension
    y = sum(x,dim);
end
