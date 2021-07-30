% Run the script to see the various outputs of POP
load testimg

% add random NaN mask to the image:
mask=(rand(size(im))<0.0005);
mask=conv2(mask,ones(21),'same');
im(mask>0)=NaN; 

%%
bParams=[2 4];
qParams=1:4;
flagpeel = 1;

s=POP(im, bParams, qParams, 'sim',flagpeel); % for sim image case
sexp=POP(im, bParams, qParams, 'exp',flagpeel); % for exp image case
figure('Position',[655 170 390 390]);
imagesc([im(:,1:end/2)./max(max(im(:,1:end/2))) s.simImage(:,end/2+1:end)./max(max( s.simImage(:,end/2+1:end)))]); axis square
title('raw vs onion peeled');

figure('Position',[1 170 650 650]);
subplot(3,3,1);imagesc(im); title('raw image');xlabel('x-pixels');ylabel('y-pixels');
subplot(3,3,2);imagesc(s.iraraw); title('iraraw');xlabel('folded angle (PPR)');ylabel('radius');
subplot(3,3,3);imagesc([0 pi/2],[0 size(s.sqp,1)],s.sqp); title('sqp');xlabel('angle [rad]');ylabel('radius');

subplot(3,3,4);imagesc(s.simImage); title('sim Image');xlabel('x-pixels');ylabel('y-pixels');
subplot(3,3,5);imagesc(s.iradecon.*(s.iradecon>0));title('iradecon');xlabel('folded angle (PPR)');ylabel('radius');
subplot(3,3,6);imagesc([0 pi/2],[0 size(s.sqpd,1)],s.sqpd.*(s.sqpd>0)); title('sqpd');xlabel('angle');ylabel('radius');

subplot(3,3,7);imagesc(sexp.expImage); title('expImage');xlabel('x-pixels');ylabel('y-pixels');
subplot(3,3,8);imagesc(sexp.iradeconExp.*(sexp.iradeconExp>0));title('iradeconExp');xlabel('folded angle (PPR)');ylabel('radius');
subplot(3,3,9);imagesc([0 pi/2],[0 size(sexp.sqpde,1)],sexp.sqpde.*(sexp.sqpde>0)); title('sqpde');xlabel('angle');ylabel('radius');


sb=size(s.Betas,1);
r=1:numel(s.PESId);
figure('Position',[654  640 785 180]);
subplot(1,sb+1,1); plot(r,s.PESId);title('PESId');xlabel('radius');ylabel('intensity');
for nsb=1:sb
      subplot(1,sb+1,1+nsb);  plot(r,s.Betas(nsb,:).*(s.Betas(1,:)>5) ); title(['\beta_{' num2str(nsb*2-2) '}']);xlabel('radius');ylabel('intensity');  
end
%% go 4D but resample because otherwise it's slow
clear be xq s3 s2 s1
[br, xq]=resample(s.Betas',1:size(s.Betas,2),64/size(s.Betas,2),5,5);
 br=br';
figure('Position',[1050 170 390 390]);
xx=[-fliplr(xq) xq(2:end)];
[X Y Z]=meshgrid(xx);
[s1 s2 s3]=beta2cart(br);

s3(s3<0.05*max(s3(:)))=nan; % take out small values for visualization
h = slice(X,Y,Z,s3, [], [], xx);
set(h, 'EdgeColor','none', 'FaceColor','interp');
% set view
axis tight; daspect([1 1 1])
view([20,25])

% normalize s3 for alpha (transparency) map
ns3=(s3(:)-min(s3(:)))./(max(s3(:))-min(s3(:)));
ns3=reshape(ns3,size(s3));
% set alpha map
for n=1:numel(h)
    set(h(n),'AlphaData',ns3(:,:,n),'FaceAlpha','interp')
end
title('s(3)');