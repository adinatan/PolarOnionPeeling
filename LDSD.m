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