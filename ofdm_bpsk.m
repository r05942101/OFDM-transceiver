function ber = ofdm_bpsk(Eb_N0)
load LTI
% define parameter
N0 = 2;
Nc = 128;
L = 30;
num = 1000;
total = num*Nc;

s_bitstream = randi([0,1],1,total);

p_bitstream = zeros(Nc,num); 
% serial to parallel
for i =1:num
    for j =1:Nc
        p_bitstream(j,i) = s_bitstream(Nc*(i-1)+j);
    end
end

Eb_N0 = 10;
U = zeros(1,Nc);
u = zeros(1,Nc);
d = sqrt(2*10.^(Eb_N0/10));
x = zeros(1,Nc+L-1);
Y = zeros(Nc,num);
y = zeros(Nc,num);
Sr = zeros(Nc,num);
Si = zeros(Nc,num);
power = zeros(Nc,num);
N1 = normrnd(0,N0/2,[Nc,num]);
N2 = normrnd(0,N0/2,[Nc,num]);
n = length(h);

for i =1:num
% bpsk modulation
    for j =1:Nc
       if (p_bitstream(j,i)==0)
           U(j) = -1*d;
       else
           U(j) = d;
       end
    end
% IDFT 
    for j =1:Nc
        temp = 0;
        for k =1:Nc
            temp = temp + U(k)*exp(2*pi*(j-1)*(k-1)*1i/128);
        end
        u(j) = temp/sqrt(Nc);
        power(j,i)= abs(u(j));
       
    end
% parallel to series & insert CP
    x(1,1:L-1) = u(1,Nc-L+2:Nc);
    x(1,L:end) = u(1,1:Nc);
 

% channel filter
    for j=1:Nc
        temp = 0;
        for k=1:n
            temp = temp + x(j+k-1)*h(n-k+1);
        end
        Y(j,i) = temp;
    end
end

% add noise
Y = Y+complex(N1,N2);


DFT_Y = zeros(Nc,num);
for i =1:num
    for j=1:Nc
        temp = 0;
        for k =1:Nc
            temp = temp+Y(k,i)*exp(-1*1i*2*pi*(j-1)*(k-1)/Nc);
        end
        DFT_Y(j,i) = temp/sqrt(Nc);
    end
end
        
H = zeros(1,Nc);
for i =1:Nc
    temp = 0;
    for j=1:n
        temp = temp+h(j)*exp(-1*1i*2*pi*(i-1)*(j-1)/Nc);
    end
    H(i) = temp;
end

rx_bitstream = zeros(1,total);
for i =1:num
    for j=1:Nc
        d1 = H(j)*d;
        d2 = H(j)*-1*d;
        if (abs(DFT_Y(j,i)-d1)<abs(DFT_Y(j,i)-d2))
            rx_bitstream((i-1)*Nc+j) = 1;
        else
            rx_bitstream((i-1)*Nc+j) = 0;
        end
    end
end

% calculate BER
error =0;
for i=1:total
    if(rx_bitstream(i)~=s_bitstream(i))
        error = error+1;
    end
end
ber = error/total
% calculate PAPR
Max = max(max(power));
Mean = mean(mean(power));
PAPR = Max/Mean
end
        
        
        

