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