%Autoregressive Gibbs Sampler%
%EFKHunt 2011%
function ARGibbsV2
l = 200; %number of missing samples
m = 500; %onset of missing data
p= 80; % Order of AR Model

rstrm = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(rstrm);
warning off

[u,x,t,N] = GenerateSignal(l,m);

figure(1); clf(1); plot(t,x); xlabel('t'); ylabel('x'); title('Restored Low-noise Sine-wave (3 iterations)'); hold on;

w = x';
z= zeros(l,1);
theta = ones(p,1); a =[0;0;0;0]; sigma = 1; %initialise the starting values

q = [l m N p]; %Vector containing all the constants which are required by the matrix generating function

for i = 1:10
%Sample z
[L,e,B,D,Y,z,wr] = Generate_Matrices(q, w, z, theta);
C_z =sigma^2*inv(D);
mu_z = -D^-1 * B' *Y; 
z = mvnrnd(mu_z,C_z);
z= z';

%Sample theta
[L,e,B,D,Y,z,wr] = Generate_Matrices(q, w,z, theta);
C_t = sigma^2*(inv(L'*L));  mu_t = (L'*L)^(-1)*L'*wr; theta  = mvnrnd(mu_t,C_t);
theta = theta';

%Sample sigma
[L,e,B,D,Y,z,wr] = Generate_Matrices(q, w, z, theta);
u = gamrnd((N-1)/2,1); sigma = (e'*e/2)^0.5/(sqrt(u));

end
z;
w(m:m+l-1)=z;
plot(w,'k-');
figure(2); plot(e,'k-');xlabel('sample number'); ylabel('Amplitude'); title('Excitation sequence'); 
sigma = sigma 

end


function [u,x,t,N] =  GenerateSignal(l,m)
s = 44100;
T = 0.0227;
N = round(s*T);
%N = 100;

f = 275.3; %Sound freq in hertz

t = [1:1:N]; %time
n = random('norm',0,0.02,1,N); %the noise signal
u = sin(f*2*pi*t/s) + n; %The noisy signal (with freq in samples s^-1

ind = linspace(m,m+l-1,l); x= u+n; x(ind) = 0; %sample x from u, with a missing part between m, m+l
end




function [L,e,B,D,Y,z,wr] = Generate_Matrices(q, w, z, theta)
%q is a vector containing l,m,N,p

l = q(1) 
m = q(2) 
N=q(3); p=q(4);

Y1= w(1:m-1); sizeY1 = size(Y1);
w(m:m+l-1) = z;
sizez =size(z);
Y2 =w(m+l:N); sizeY2 = size(Y2);
sizew = size(w);
%'De-augment' the w vector into its component parts

Y = [Y1;Y2];
%K_diag = zeros(N,p);


L_diag = zeros(size(w),size(w)+p-1);

L = zeros(size(w),size(p));
for i=1:size(w)-1;
    L_diag(:,i+p-1) = w(i);
end
L_diag;
L = zeros(size(w)-1,p);
L_diag_ind = [p-1:-1:-sizew+2];
L = spdiags(L_diag,L_diag_ind,L);
%L(size(L)-p-2: size(L),:)  = [];
%L = L
L=full(L)

r = zeros(1,size(w)); t2 = zeros(1,size(w)-size(theta)-1); K = toeplitz([1 -theta' t2],r);

KTK = K'*K;
B1 = KTK(1:sizeY1,sizeY1+1:sizeY1+sizez);
B2 = KTK(sizeY1+sizez+1:sizeY1+sizez+sizeY2,sizeY1+1:sizeY1+sizez); %Extract B2
D = KTK(sizeY1+1:sizeY1+sizez,sizeY1+1:sizeY1+sizez);
B = [B1;B2];

sizew = size(w); sizeL = size(L); sizet = size(theta);
wr = w; wr(1)= [];
e = wr - L*theta;


end
