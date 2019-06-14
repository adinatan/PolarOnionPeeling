function ff=bld_lut(r,fsize,N)
% The function builds the basis set of radial delta functions used in POP.m
% The code simulates a 3D sphere of radius r with uniformly distributed 
% random  number of events. Then project the 3D sphere into a Cartesian 
% 2D image s(r), similar to what in measured in experiments. Then we obtain
% the radial projection of s(r) using the same triangular polar 
% representation that is used in the POP code. I shared this code for the
% users reference and for it's pedagogical value. The user is warned that
% creating the full lut.mat file with the appropriate accuracy (N=1e8) for
% images of size 1024x1024 (i.e. 1<=r<=512) takes several hours. That is
% the reason lut.mat is added  as a separate file.
%
% Inputs:
%  r      - The radius of the sphere that is projected to 2D.
%  fsize  - The final size of the basis set needed.
%  N      - # of simulated events per radius. In the lut.mat file N=1e8. 
% 
%
% Outputs:
%  ff     -  a vector of the radial projection 
%
%   Example:
%       ff=bld_lut(400,1025,1e7);
%       plot(ff);
%
%   Comments \ improvements are welcomed
%   Adi Natan (natan@stanford.edu)
 
 %r=r+1;
x=@(r,theta,phi) r.*sin(theta).*sin(phi);  
z=@(r,theta) r.*cos(theta);
 
t=rand(N,1)*pi;   % theta
p=rand(N,1)*2*pi; % phi
 
% 3d weighted histogram as function of r
s=@(r) hist3w( x(r,t,p),  z(r,t), sin(t), -r-0.5:r+0.5, -r-0.5:r+0.5 );
%s=@(r) hist3w( x(r,t,p),  z(r,t), sin(t), -r-1:r+1, -r-1:r+1 );

% get radial intensity and normalize
ff=radialI(s(r)+eps);
ff(fsize)=0;
ff=ff./max(ff);

function B=radialI(im)
% Cartesian to polar coordinates with quadrant symmetry
%
% Inputs:
%  im       - a  map image, the image is assumed to be square
%             matrix, centered and corrected for deformations, tilt, etc.
%             the image may contain NaN, and the function will ignore them.

qParams=1:4; % include all image unless stated otherwise.
x0=ceil(size(im,1)/2); y0=ceil(size(im,2)/2); % Assuming image is centered
RR=(0:size(im)/2);
PPR=(floor(0.5*pi*(RR+1))-1); % calc the  # of pixels per radius
AngleInc = (0.5*pi./PPR'); % angle increment per radius
AngleInc(1)=0; % avoid inf at origin

%set radius of square matrix around center
L = min([x0,y0,(size(im) - [x0,y0])]);

 
%create  image quadrants
if mod(size(im,1),2)==1 %case for even image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0:y0+L-1);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0:x0+L-1   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0:x0+L-1   ,y0:y0+L-1);
else %case for odd image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0+1:y0+L);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0+1:x0+L   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0+1:x0+L   ,y0+1:y0+L);
end

%add image quadrants.  NaN's or zeros are not counted if there are values
%from other quadrant, multiple values are then averaged. This assumes that
%each quadrant has the same response or sensitivity, if not this need to be
%corrected
 
a4=nansum(Q(:,:,qParams),3)./sum( ~isnan(Q(:,:,qParams)).*(Q(:,:,qParams)~=0),3);
 
    
ira = zeros(L,PPR(L)); % initialize the  matrix
ira(1,1) = a4(1,1);      % origin pixel remains the same
B(L)=0;

%% creating the 2d triangular array polar image 
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
    % gather intensity per radius per angle (ira)
    ira(r,1:npr+1) = xd.*yd.*a4(xc+(yc-1)*L);
    ira(r,2:npr) = ira(r,2:npr) + ...
        xd(2:npr).*(1-yd(2:npr)).*a4(xc(2:npr)+(yc(2:npr)-1+(-1).^(yp(2:npr)<yc(2:npr)))*L) + ...
        (1-xd(2:npr)).*yd(2:npr).*a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+yc(2:npr)*L) + ...
        (1-xd(2:npr)).*(1-yd(2:npr)).*a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+L*(yc(2:npr)+(-1).^(yp(2:npr)<yc(2:npr))));
    
   
   %y = ira(r,qp+1)-eps; 
   B(r)=sum(ira(r,qp+1));%/sum(~isnan(y))*(npr+1); 
end


function f  = hist3w(x,y,weights,xedges,yedges)
%hist3 with weights,  all inputs should be vectors
[~,xbin] = histc(x,xedges);
[~,ybin] = histc(y,yedges);
 
%xbin, ybin zero for out of range values
% (see the help of histc) force this event to the
% first bins
xbin(xbin==0) = inf;
ybin(ybin==0) = inf;
xnbin = numel(xedges);
ynbin = numel(yedges);
 
xy = xbin * ynbin + ybin;
 
[xyuni,id] = unique(xy);
xyuni(end) = []; % remove Inf bin
id(end) = [];
hstres = histc(xy,xyuni);
hstres = hstres.*weights(id); % weigh the histogram
f(ynbin,xnbin)=0; % preallocate memory
f(xyuni-ynbin) = hstres;
