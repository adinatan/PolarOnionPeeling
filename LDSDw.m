function [betas varargout]=LDSDw(vec, bParams, w)
warning('off','all')

% Apply Legendre decomposition for a single data vector spanning [0,2*pi]
% using weigted least squares (WLS)
%
% Inputs:
%  vec      - a vector that contain information that spans 0 to 2*pi
%  bParams  - a vector with the (even) beta parameter (legenadre order)
%             to be used in the fitting procedure (excluding 0th-order
%             which is always included), e.g. [2 4] for beta_2 & beta_4.
%  w        - a vector of error weights per element in vec, usualy the by variance 
%
% Outputs:
%  betas      - a vector containing each beta_n such that positive
%               orders are normalized by beta_0 (intensity of vec),
%               for example: betas(1) is beta_0, betas(2) is beta_2/beta_0,
%               betas(3) is beta4/beta_0, ...
%  ese        - a vector containing the estimated standard error of the WLS
%  reco_vec   - the reconstructed vector using the betas found
%
%   Ver 1 (2019-11-14)
%   Adi Natan (natan@stanford.edu)
%

% % example:
% % generate random vec up to beta 4 + noise

% N=64; % # of bins in angle from 0 to 2*pi
% rv=rand(1,3); % generate random coeffiients 
% w=randi(1e3,1,N); % random weights vector
% werr=1./(sqrt(w));
% 
% vec_noisy=100*(rv(1)-0.5)+... % baseline
%     (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
%     (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
%     werr.*randn(1,N); % noise term
% 
% vec_gt=100*(rv(1)-0.5)+... % baseline
%     (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
%     (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
%     0*randn(1,N); % no noise term
% 
% 
% %apply Legendre decomposition up to a higher order and see that indeed
% %the higher order doesnt contribute
% [betas, ese, reco_vec]=LDSDw(vec_noisy, 2:2:6, w);
% [betas_gt]=LDSDw(vec_gt, 2:2:6);
% 
% % plot it all
% errorbar(linspace(0,2*pi,N),vec_noisy,werr,'o','Color',[0.5 0.5 0.5]); hold on;
% plot(linspace(0,2*pi,N),vec_gt,'k')
% plot(linspace(0,2*pi,N),reco_vec,'r');
% legend('weighted noisy vec input','ground truth','WLS fit');
%   
% St{1}= ['\bf\color{red} ', '\beta_n from ground truth:', ' \rm\color{black}'];
% Sc{1}= ['\bf\color{red} ', '\beta_n from WLS fit:', ' \rm\color{black}'];
% for n=1:numel(betas)
%     if n>1
%         St{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str(betas_gt(n)*betas_gt(1)) ];
%         Sc{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str(betas(n)*betas(1)) '\pm'  num2str(ese(n))];
%     else
%         St{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str( betas_gt(n)) ];
%         Sc{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str( betas(n)) '\pm'  num2str(ese(n))];
%     end
% end
% St{n+1}='';
% h = gca; % get the current axis;
% h.Position(3) = 0.5;
% annotation('textbox', [0.65, 0.825, 0.1, 0.1], 'String', cat(2,St,Sc))

%% defaults
if (nargin < 2);  bParams=[2 4]; w=ones(size(vec)); end
if (nargin < 3);                 w=ones(size(vec)); end

% Check that the beta Parameters are in the scope of the code
if any(mod(bParams,2)) || any(bParams<=0)
    error('Only even positive beta parameters supported!');
end


% helper function to replace the use the leg in the original code:
paren = @(x, varargin) x(varargin{:});



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
notnan=~isnan(y);
B(1)=nansum(y)/sum(notnan)*(npr+1);

% one fit coefficient for each B param
fitCoefs = ones(numel(bParams) + 1, sum(notnan));
for ii=(1:numel(bParams))
    % assign relevant Legendre polynomials to fitCoef
    fitCoefs(ii+1,:) = paren( legendre(bParams(ii), cos(AngleInc*qp(notnan))), 1,1:sum(notnan));
    B(ii+1) = y(notnan)*fitCoefs(ii+1,:)'; % fill B vector
end

[bw,ese] = lscov(fitCoefs',y(notnan)',w(notnan)); % weighted LS
 
 
betas(1)=bw(1);
betas(2:end)=bw(2:end)./bw(1);


%% reconstruct signal
if nargout>1
  
    alphaMat = (AngleInc * (0:PPR));
    rrpCosAlphaMat = ones(1,PPR+1).*cos(alphaMat);
    
    clear beta_contrw
    for ii=1:numel(bParams) %  add each beta contribution
        beta_contrw(:,ii)= bw(ii+1)./bw(1)*paren( legendre(bParams(ii), rrpCosAlphaMat), 1,1:numel(y));
    end
    
    factMatw = bw(1).*(ones(1, PPR+1) + sum(beta_contrw,2));
    reco_vecw(1:PPR+1)=factMatw(1:PPR+1,1);% /sqrt(r);
    
    varargout{1} = ese; %estimated_standard_error
    varargout{2}= reco_vecw; %reconstructed vector
end
