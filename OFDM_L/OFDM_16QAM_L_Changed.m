function BER=OFDM_16QAM_L_Changed(L)
%EbN0=0
load LTI;
G=1000;
N=128*G*4;
Nc = 128;
%Generate random bits (1 or 0)
r = randi(0:1,N,1);
ra = zeros(128*4,1,G);
for w=1:G
    ra(:,1,w)=r((w-1)*128*4+1:(w)*128*4);
end



%% Converting bits to symbol under 8PSK modulation
[row,col]=size(ra);
EbN0 = 10;
d = sqrt(0.8*10.^(0.1*EbN0));
R = zeros(128,1,G);
S = R;
X = zeros(Nc+L-1,1,G);
Y = zeros(128,1,G);
PAPR=zeros(128,1,G);
for w=1:G
    for p=1:4:(row-3)
    %0 1 3 4 6 7 5 4
        if(ra(p,1,w)==0 && ra(p+1,1,w)==0 && ra(p+2,1,w)==0  && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(-3*d,-3*d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==0 && ra(p+2,1,w)==0 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(-3*d,-d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==0 && ra(p+2,1,w)==1 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(-3*d,d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==0 && ra(p+2,1,w)==1 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(-3*d,3*d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==1 && ra(p+2,1,w)==0 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(-d,-3*d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==1 && ra(p+2,1,w)==0 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(-d,-d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==1 && ra(p+2,1,w)==1 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(-d,d);
        elseif (ra(p,1,w)==0 && ra(p+1,1,w)==1 && ra(p+2,1,w)==1 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(-d,3*d);
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==1 && ra(p+2,1,w)==0 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(d,-3*d);   
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==1 && ra(p+2,1,w)==0 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(d,-d);        
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==1 && ra(p+2,1,w)==1 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(d,d);        
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==1 && ra(p+2,1,w)==1 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(d,3*d);        
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==0 && ra(p+2,1,w)==0 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(3*d,-3*d);   
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==0 && ra(p+2,1,w)==0 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(3*d,-d); 
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==0 && ra(p+2,1,w)==1 && ra(p+3,1,w)==1)
            R((p-1)/4+1,1,w)=complex(3*d,d);      
        elseif (ra(p,1,w)==1 && ra(p+1,1,w)==0 && ra(p+2,1,w)==1 && ra(p+3,1,w)==0)
            R((p-1)/4+1,1,w)=complex(3*d,3*d);
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
    
    %% Cyclic Prefix Insertion

   X(1:(L-1),1,w)=S((Nc-L+2):Nc,1,w);
   X(L:end,1,w) = S(1:Nc,1,w);
   PAPR(:,1,w)= abs(S(:,1,w));
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

        a1=sqrt(128)*H(b,1)*complex(-3*d,-3*d);
        a2=sqrt(128)*H(b,1)*complex(-3*d,-d);
        a3=sqrt(128)*H(b,1)*complex(-3*d,d);
        a4=sqrt(128)*H(b,1)*complex(-3*d,3*d);
        a5=sqrt(128)*H(b,1)*complex(-d,-3*d);
        a6=sqrt(128)*H(b,1)*complex(-d,-d);
        a7=sqrt(128)*H(b,1)*complex(-d,d);
        a8=sqrt(128)*H(b,1)*complex(-d,3*d);
        a9=sqrt(128)*H(b,1)*complex(d,-3*d);
        a10=sqrt(128)*H(b,1)*complex(d,-d);
        a11=sqrt(128)*H(b,1)*complex(d,d);
        a12=sqrt(128)*H(b,1)*complex(d,3*d);
        a13=sqrt(128)*H(b,1)*complex(3*d,-3*d);
        a14=sqrt(128)*H(b,1)*complex(3*d,-d);
        a15=sqrt(128)*H(b,1)*complex(3*d,d);
        a16=sqrt(128)*H(b,1)*complex(3*d,3*d);

        d1=abs(V(b,1,a)-a1);
        d2=abs(V(b,1,a)-a2);
        d3=abs(V(b,1,a)-a3);
        d4=abs(V(b,1,a)-a4);
        d5=abs(V(b,1,a)-a5);
        d6=abs(V(b,1,a)-a6);        
        d7=abs(V(b,1,a)-a7);
        d8=abs(V(b,1,a)-a8);
        d9=abs(V(b,1,a)-a9);
        d10=abs(V(b,1,a)-a10);
        d11=abs(V(b,1,a)-a11);
        d12=abs(V(b,1,a)-a12);
        d13=abs(V(b,1,a)-a13);
        d14=abs(V(b,1,a)-a14);
        d15=abs(V(b,1,a)-a15);
        d16=abs(V(b,1,a)-a16);
        

        D1=[d1 d2 d3 d4 d5 d6 d7 d8];
        D2=[d9 d10 d11 d12 d13 d14 d15 d16];
        mm=min(D1);
        nn=min(D2);
        zz=min(mm,nn);
        switch(zz)
            case d1,
                T(b,1,a)=0;
            case d2,
                T(b,1,a)=1;
            case d3,
                T(b,1,a)=3;
            case d4,
                T(b,1,a)=2;
            case d5,
                T(b,1,a)=4;
            case d6,
                T(b,1,a)=5;
            case d7,
                T(b,1,a)=7;
            case d8,
                T(b,1,a)=6;
            case d9,
                T(b,1,a)=12;
            case d10,
                T(b,1,a)=13;
            case d11,
                T(b,1,a)=15;
            case d12,
                T(b,1,a)=14;
            case d13,
                T(b,1,a)=8;
            case d14,
                T(b,1,a)=9;
            case d15,
                T(b,1,a)=11;
            otherwise,
                T(b,1,a)=10;
        end
        
    end
end

%% Converting symbol to binary bits
O=zeros(128,4,G);
O2=zeros(128*4,1,G);
for u=1:G
    O(:,:,u)=de2bi(T(:,:,u),'left-msb');
end

for u=1:G
    for p=1:128
        for q=1:4
        O2((p-1)*4+q,1,u)=O(p,q,u);
        end
    end
end


%% Calculate BER
error = 0;
for u=1:G
    for p=1:128*4
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