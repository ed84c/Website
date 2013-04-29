
%Expectation Maximisation Algorithm%
%EFKHunt 2011%
function EM
l = 200; %number of missing samples
m = 500; %onset of missing data
p= 80; % Order of AR Model

rstrm = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(rstrm);
warning off

[u,x,t,N] = GenerateSignal(l,m);

figure(1); clf(1); plot(t,x); xlabel('t'); ylabel('x'); title('Restored Quiet Sine-wave'); axis([0 1000 -1.1 1.1]); hold on;

w = x';

theta = ones(p,1); a =[0;0;0;0]; sigma = 1; %initialise the starting values

q = [l m N p]; %Vector containing all the constants which are required by the matrix generating function

z= zeros(l,1);

for i = 1:14
    
[L,e,B,D,Y,wr,Y1,Y2, T] = Generate_Matrices(q, w, z, theta);

v = wr; M = L; size(M);
size(v); phi_peak = (M'*M)^-1 * M' * v;
theta_peak = (L'*L)^-1 * L' * wr; size(phi_peak);

[L,e,B,D,Y,wr,Y1,Y2, T] = Generate_Matrices(q, w, z, phi_peak);
sigma = sqrt(e'*e/1000);

qp = conv2(T(:,20),[Y1;zeros(l,1);Y2]);
qp = qp(m:m+l-1);

size(B'); size(Y); size(q); size(B'*Y); size(qp);
z = -(sigma^2*qp+B'*Y);
size(T); size(D); size(z);
z = ((sigma^2*T+D)^-1) * z ;

end 


w(m:m+l-1)=z;
plot(w,'k-');
figure(2); plot(e,'k-');xlabel('sample number'); ylabel('Amplitude'); title('Excitation sequence'); 
sigma = sigma;

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




function [L,e,B,D,Y,wr,Y1,Y2,T] = Generate_Matrices(q, w, z, theta)
%q is a vector containing l,m,N,p

l = q(1); 
m = q(2);
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
L=full(L);

r = zeros(1,size(w)); t2 = zeros(1,size(w)-size(theta)-1); K = toeplitz([1 -theta' t2],r);

KTK = K'*K;
B1 = KTK(1:sizeY1,sizeY1+1:sizeY1+sizez);
B2 = KTK(sizeY1+sizez+1:sizeY1+sizez+sizeY2,sizeY1+1:sizeY1+sizez); %Extract B2
D = KTK(sizeY1+1:sizeY1+sizez,sizeY1+1:sizeY1+sizez);
B = [B1;B2];

sizew = size(w); sizeL = size(L); sizet = size(theta);
wr = w; wr(1)= [];
e = wr - L*theta;

Tl = (L'*L)^(-1);
len_z = size(z,1);
T = zeros(len_z);
size(T);

j = 1;
linv = (L'*L)^(-1);
[ml,nl] = size((L'*L)^(-1));
BT = zeros(len_z,(ml-1)*2);

i1 = -(nl-1);
in = ml-1;
for i= i1:in
    Td = diag(Tl,i);
    Td = sum(Td);
    BT(1:len_z,j) = Td;
    j = j + 1;
end
BT;
dT = linspace(-(nl-1),ml-1,2*(ml-1)+1); 
T = spdiags(BT,dT,T); T = full(T);
%Generates the T Matrix
end
