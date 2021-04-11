function BER=OFDM_BPSK_L_Changed(L)
%EbN0=6
load LTI;
G=1000;
N=128*G;
Nc = 128;

%Generate random bits (1 or 0)
r = randi(0:1,N,1);
ra = zeros(128,1,G);
for w=1:G
    ra(:,1,w)=r((w-1)*128+1:(w)*128);
end
EbN0 = 10;
d = sqrt(2*10.^(0.1*EbN0));
R = zeros(128,1,G);
Sr = R;
Si = R;
Tr = Sr;
Ti = Si;
X = zeros(Nc+L-1,1,G);
Y1 = zeros(Nc+L-1,1,G);
Y = zeros(128,1,G);

N0=2;
N1=normrnd(0,N0/2,[Nc,1,G]);
N2=normrnd(0,N0/2,[Nc,1,G]);
PAPR=zeros(128,1,G);
for w=1:G
    %% BPSK Modulation
    for j=1:128
        if (ra(j,1,w) == 0)
            R(j,1,w) = -1*d;
        else
            R(j,1,w) = d;
        end
    end
    
    %% IDFT
    
    for j=0:127
        sr = 0;
        si = 0;
        for k=0:127
            sr = sr + R((k+1),1,w)*cos(2*pi*j*k/128);
            si = si + R((k+1),1,w)*sin(2*pi*j*k/128);
        end
        Sr(j+1,1,w) = 1/sqrt(128)*sr;
        Si(j+1,1,w) = 1/sqrt(128)*si;
        PAPR(j+1,1,w)= (Sr(j+1,1,w)^2 + Si(j+1,1,w)^2);
    end 
    %% Cyclic Prefix Insertion
   xr = zeros(Nc+L-1,1,G);
   xi = zeros(Nc+L-1,1,G);
   
   xr(1:(L-1),1,w)=Sr((Nc-L+2):Nc,1,w);
   xi(1:(L-1),1,w)=Si((Nc-L+2):Nc,1,w); 
   
   xr(L:end,1,w) = Sr(1:Nc,1,w);
   xi(L:end,1,w) = Si(1:Nc,1,w);
   X(:,:,w) = complex(xr(:,:,w),xi(:,:,w));
    %% LTI Filter
    for m=1:Nc+L-1
        tmp=0;
        for n=1:length(h)
            if(m>=n)
            tmp = tmp + X(m-n+1,1,w)*h(n);
            end
        end
        Y1(m,1,w) = tmp;
    end    
end
    %% CP Removal
    for p=L:Nc+L-1
        Y(p-L+1,1,:)=Y1(p,1,:);
    end
%% Add noise
Y=Y+complex(N1,N2);

%% DFT
V = Y;
for a = 1:G
    for k = 1:Nc
        tmp = 0;
        for n=0:Nc-1
            dft = exp(-2*pi*i*(k-1)*n/128);
            tmp = tmp+Y(n+1,1,a)*dft;
        end
        tmp=tmp/sqrt(128);
        V(k,1,a) = tmp;
    end
end

%% Detection
H = zeros(128,1);
hh = zeros(128,1);
hh(1:30) = h(:);
%H(1:128,1) = hh(:);

for m = 1:Nc
    tmp = 0;
    for n = 0:(128-1)
        dd = exp(-2*pi*i*(m-1)*n/128);
        tmp = tmp+hh(n+1)*dd;
    end
    tmp = tmp/sqrt(128);
    H(m,1)=tmp;
end

T = zeros(Nc,1,G);

for a=1:G
    for b=1:Nc
        d1=sqrt(128)*H(b,1)*d;
        d2=sqrt(128)*H(b,1)*(-d);
        if(abs(V(b,1,a)-d1) <= abs(V(b,1,a)-d2))
            T(b,1,a) = 1;
        else
            T(b,1,a) = 0;
        end
    end
end


%% Calculate BER
error = 0;
for a=1:G
    for b=1:Nc
        if(T(b,1,a) ~= ra(b,1,a))
            error = error + 1;
        end
    end
end
BER = error / N;
end



