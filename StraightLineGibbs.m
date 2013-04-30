%Straight Line Gibbs Sampler%
%EFKHunt 2011%
function StraightLineGibbs
l = 30; %number of missing samples
m = 40; %onset of missing data
N = 100; % length of data


%rstrm = RandStream('mt19937ar','Seed',2)
%RandStream.setGlobalStream(rstrm)

[u,x,t] = GenerateSignal(l,m,N);


figure(1); clf(1); plot(t,x); xlabel('t'); ylabel('u'); title('Noisy straight line signal with missing data. 40 iterations.'); hold on;
d = u; %sets the sampled data vector equal to the data vector generated initialy
D1 = sum(d); D2 = sum(d.*d); Xd = sum(t.*d); X1 = sum(t); X2 = sum(t.*t); %Sets up the constants required for the sampling

mu=20; sigma=2; c =3; m=4; %initialise the starting values

for i = 1:20
D1 = sum(d); D2 = sum(d.*d); Xd = sum(t.*d); X1 = sum(t); X2 = sum(t.*t); 
mu_m = (Xd-c*X1)/X2; sigma_m = sigma/sqrt(X2); m = normrnd(mu_m,sigma_m); %iteration on finding m 

D1 = sum(d); D2 = sum(d.*d); Xd = sum(t.*d); X1 = sum(t); X2 = sum(t.*t); 
mu_c = (D1-m*X1)/N; sigma_c = sigma/sqrt(N); c = normrnd(mu_c,sigma_c); %iteration finding c
u-m*t-c;

alpha = N/2; beta = 1/2*sum((u-m*t-c).^2)
c_g = 2*beta^alpha /gamma(alpha)
sigma = SRIG(alpha,beta,c_g);
end
m =m
c_f =c 
sigma_f =sigma
plot(t,t*m+c + random('norm',0,sigma,1,N),'k-');




end


function [u, x,t] =  GenerateSignal(l,m,N)
t = [0:1:N-1]; %time
n = random('norm',0,1,1,N); %the noise signal
u = 0.5*t+n+1; %The noisy signal
%figure(1); clf(1); plot(t,u); xlabel('t'); ylabel('u'); title('Noisy sine signal');

ind = linspace(m,m+l,l+1); x= u; x(ind) = NaN; %sample x from u, with a missing part between m, m+l
end

function theta= SRIG(alpha,beta,c)
c = 2*beta^alpha/gamma(alpha);
x =[0:0.01:10];
p_SRIG = @(alpha, beta, c, x) c*x.^(-2*alpha+1).*exp(-beta./x.^2); %this is the pdf for the Square Root Inverted Gamma (SRIG) distribution


q = normpdf(x,0,1);
k = 0.01;


%figure(2); plot(x,p_SRIG(alpha, beta, c, x)); hold on; 
while sum(k*q < p_SRIG(alpha, beta, c, x)) > 0
    k = k + 0.1;
%figure(2); plot(x,k*q,'k-'); title('Square Root Inverted gamma PDF, and c*Guassian pdf for increasing c.')      
end 
%Find a k such that kq is always greater than p_SRIV ie. a k for  k * normal dist.
%pdf where k * normal dist 'envelopes' over the top of the SRIG. The sum
%is used in the test since '>' returns a vector for each element of true/
%false for greater than.

reject = true;
r=0;
while reject == true;
r = r +1;
    theta = normrnd(0,1); %sample theta from the normal distribution
    u = unifrnd(0,normpdf(theta,0,1)*k); %sample a uniform deviate from the inteval from [0,kq(theta)]
    p_SRIG(alpha, beta, c, theta);
    if u < p_SRIG(alpha, beta, c, theta);
        reject = false;
    end 
    %Reject this theta if u> p(theta);
    
    if r>100
        error('Stuck in while loop')
    end
end 

end