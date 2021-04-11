function BER=OFDM_8PSK_L_Changed(L)
%EbN0=0
load LTI;
G=1000;
N=128*G*3;
Nc = 128;

%Generate random bits (1 or 0)
r = randi(0:1,N,1);
ra = zeros(128*3,1,G);
for w=1:G
    ra(:,1,w)=r((w-1)*128*3+1:(w)*128*3);
end



%% Converting bits to symbol under 8PSK modulation
[row,col]=size(ra);
EbN0 = 10;
d = sqrt(6*10.^(0.1*EbN0));
R = zeros(128,1,G);
S = R;
X = zeros(Nc+L-1,1,G);
Y = zeros(128,1,G);
PAPR=zeros(128,1,G);
for w=1:G
    for e=1:3:(row-2)
        if(ra(e,1,w)==0 && ra(e+1,1,w)==0 && ra(e+2,1,w)==0)
        R((e-1)/3+1,1,w)=complex(d,0);
        elseif (ra(e,1,w)==0 && ra(e+1,1,w)==0 && ra(e+2,1,w)==1)
        R((e-1)/3+1,1,w)=complex(d/sqrt(2),d/sqrt(2));
        elseif (ra(e,1,w)==0 && ra(e+1,1,w)==1 && ra(e+2,1,w)==1)
        R((e-1)/3+1,1,w)=complex(0,d);
        elseif (ra(e,1,w)==0 && ra(e+1,1,w)==1 && ra(e+2,1,w)==0)
        R((e-1)/3+1,1,w)=complex(-d/sqrt(2),d/sqrt(2));
        elseif (ra(e,1,w)==1 && ra(e+1,1,w)==1 && ra(e+2,1,w)==0)    
        R((e-1)/3+1,1,w)=complex(-d,0);
        elseif (ra(e,1,w)==1 && ra(e+1,1,w)==1 && ra(e+2,1,w)==1)
        R((e-1)/3+1,1,w)=complex(-d/sqrt(2),-d/sqrt(2)); 
        elseif (ra(e,1,w)==1 && ra(e+1,1,w)==0 && ra(e+2,1,w)==1)
        R((e-1)/3+1,1,w)=complex(0,-d);
        elseif (ra(e,1,w)==1 && ra(e+1,1,w)==0 && ra(e+2,1,w)==0)%111
        R((e-1)/3+1,1,w)=complex(d/sqrt(2),-d/sqrt(2));
        end
    end
   
    %% IDFT
    for q=0:127
        s = 0;
        for k=0:127
            s = s + R((k+1),1,w)*exp(2*pi*q*i*k/128);
        end
        S(q+1,1,w) = 1/sqrt(128)*s;
        
    end
    PAPR(:,1,w)= abs(S(:,1,w));
    %% Cyclic Prefix Insertion

   X(1:(L-1),1,w)=S((Nc-L+2):Nc,1,w);
   X(L:end,1,w) = S(1:Nc,1,w);

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
N0=2;
N1=normrnd(0,N0/2,[Nc,1,G]);
N2=normrnd(0,N0/2,[Nc,1,G]);
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

%% Filter H in Frequency Domain
H = zeros(128,1);
hh = zeros(128,1);
hh(1:30) = h(:);

for m = 1:Nc
    tmp = 0;
    for n = 0:(128-1)
        dd = exp(-2*pi*i*(m-1)*n/128);
        tmp = tmp+hh(n+1)*dd;
    end
    tmp = tmp/sqrt(128);
    H(m,1)=tmp; %H is the "h" in frequency domain
end
%% Detection
T = zeros(Nc,1,G);
for a=1:G
    for b=1:Nc

        a1=sqrt(128)*H(b,1)*complex(d,0);
        a2=sqrt(128)*H(b,1)*complex(d/sqrt(2),d/sqrt(2));
        a3=sqrt(128)*H(b,1)*complex(0,d);
        a4=sqrt(128)*H(b,1)*complex(-d/sqrt(2),d/sqrt(2));
        a5=sqrt(128)*H(b,1)*complex(-d,0);
        a6=sqrt(128)*H(b,1)*complex(-d/sqrt(2),-d/sqrt(2));
        a7=sqrt(128)*H(b,1)*complex(0,-d);
        a8=sqrt(128)*H(b,1)*complex(d/sqrt(2),-d/sqrt(2));

        d1=abs(V(b,1,a)-a1);
        d2=abs(V(b,1,a)-a2);
        d3=abs(V(b,1,a)-a3);
        d4=abs(V(b,1,a)-a4);
        d5=abs(V(b,1,a)-a5);
        d6=abs(V(b,1,a)-a6);        
        d7=abs(V(b,1,a)-a7);
        d8=abs(V(b,1,a)-a8);

        D=[d1 d2 d3 d4 d5 d6 d7 d8];
        mm=min(D);
        
        switch(mm)
            case d1,
                T(b,1,a)=0;
            case d2,
                T(b,1,a)=1;
            case d3,
                T(b,1,a)=3;
            case d4,
                T(b,1,a)=2;
            case d5,
                T(b,1,a)=6;
            case d6,
                T(b,1,a)=7;
            case d7,
                T(b,1,a)=5;
            otherwise,
                T(b,1,a)=4;
        end
        
    end
end

%% Converting symbol to binary bits
O=zeros(128,3,G);
O2=zeros(384,1,G);
for u=1:G
    O(:,:,u)=de2bi(T(:,:,u),'left-msb');
end

for u=1:G
    for p=1:128
        for q=1:3
        O2((p-1)*3+q,1,u)=O(p,q,u);
        end
    end
end


%% Calculate BER
error = 0;
for u=1:G
    for p=1:384
            if(O2(p,1,u) ~= ra(p,1,u))
                error = error +1;
            end
    end
end
  
BER = error / N;

%% Calculate PAPR
MAX1 = max(PAPR);
MEAN1 = mean(PAPR);
MAX = max(MAX1);
MEAN = mean(MEAN1);

PAPRr = MAX / MEAN;


        
end
