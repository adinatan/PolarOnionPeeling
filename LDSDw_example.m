close all
clear all
N=2^6; % # of bins in angle from 0 to 2*pi
rv=rand(1,3); % generate random coeffiients 
w=randi(3e3,1,N); % random weights vector
werr=1./(sqrt(w));

vec_noisy=100*(rv(1)-0.5)+... % baseline
    (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
    (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
    werr.*randn(1,N); % noise term

vec_gt=100*(rv(1)-0.5)+... % baseline
    (rv(2)-0.5)*cos(linspace(0,2*pi,N)).^2+... % L2
    (rv(3)-0.5)*cos(linspace(0,2*pi,N)).^4+... % L4
    0*randn(1,N); % no noise term


%apply Legendre decomposition up to a higher order and see that indeed
%the higher order doesnt contribute
[betas, ese, reco_vec]=LDSDw(vec_noisy, 2:2:6, w);
[betas_gt]=LDSDw(vec_gt, 2:2:6);

% plot it all
errorbar(linspace(0,2*pi,N),vec_noisy,werr,'o','Color',[0.5 0.5 0.5]); hold on;
plot(linspace(0,2*pi,N),vec_gt,'k')
plot(linspace(0,2*pi,N),reco_vec,'r');
legend('weighted noisy vec input','ground truth','WLS fit');
  
St{1}= ['\bf\color{red} ', '\beta_n from ground truth:', ' \rm\color{black}'];
Sc{1}= ['\bf\color{red} ', '\beta_n from WLS fit:', ' \rm\color{black}'];
for n=1:numel(betas)
    if n>1
        St{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str(betas_gt(n)*betas_gt(1)) ];
        Sc{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str(betas(n)*betas(1)) '\pm'  num2str(ese(n))];
    else
        St{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str( betas_gt(n)) ];
        Sc{n+1}=['\beta_' num2str(2*n-2) ' = ' num2str( betas(n)) '\pm'  num2str(ese(n))];
    end
end
St{n+1}='';
h = gca; % get the current axis;
h.Position(3) = 0.5;
annotation('textbox', [0.65, 0.825, 0.1, 0.1], 'String', cat(2,St,Sc))
 