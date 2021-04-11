function ber = ofdm_8psk(Eb_N0)
load LTI
% define parameter
N0 = 2;
Nc = 128;
L = 30;
num = 1000;
total = num*Nc*3;

s_bitstream = randi([0,1],1,total);

p_bitstream = zeros(Nc*3,num); 
% serial to parallel
for i =1:num
    for j =1:Nc*3
        p_bitstream(j,i) = s_bitstream(Nc*3*(i-1)+j);
    end
end


u = zeros(1,Nc);
d = sqrt(6*10.^(Eb_N0/10));
x = zeros(1,Nc+L-1);
Y = zeros(Nc,num);
y = zeros(Nc,num);
Sr = zeros(Nc,num);
Si = zeros(Nc,num);
power = zeros(Nc,num);
N1 = normrnd(0,N0/2,[Nc,num]);
N2 = normrnd(0,N0/2,[Nc,num]);
n = length(h);
row = size(p_bitstream);
for i =1:num
% bpsk modulation
    U = zeros(1,Nc);
    for j =1:3:(row-2)
        if(p_bitstream(j,i)==0 && p_bitstream(j+1,i)==0 && p_bitstream(j+2,i)==0)
            U((j-1)/3+1) = complex(d,0);
        elseif (p_bitstream(j,i)==0 && p_bitstream(j+1,i)==0 && p_bitstream(j+2,i)==1)
            U((j-1)/3+1) = complex(d/sqrt(2),d/sqrt(2));
        elseif (p_bitstream(j,i)==0 && p_bitstream(j+1,i)==1 && p_bitstream(j+2,i)==1)
            U((j-1)/3+1) = complex(0,d);
        elseif (p_bitstream(j,i)==0 && p_bitstream(j+1,i)==1 && p_bitstream(j+2,i)==0)
            U((j-1)/3+1) = complex(-d/sqrt(2),d/sqrt(2));
        elseif (p_bitstream(j,i)==1 && p_bitstream(j+1,i)==1 && p_bitstream(j+2,i)==0)   
            U((j-1)/3+1) = complex(-d,0);
        elseif (p_bitstream(j,i)==1 && p_bitstream(j+1,i)==1 && p_bitstream(j+2,i)==1) 
            U((j-1)/3+1) = complex(-d/sqrt(2),-d/sqrt(2)); 
        elseif (p_bitstream(j,i)==1 && p_bitstream(j+1,i)==0 && p_bitstream(j+2,i)==1) 
            U((j-1)/3+1) = complex(0,-d);
        elseif (p_bitstream(j,i)==1 && p_bitstream(j+1,i)==0 && p_bitstream(j+2,i)==0) 
            U((j-1)/3+1) = complex(d/sqrt(2),-d/sqrt(2));
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
    for j=1: Nc
        d1 = H(j)*complex(d,0);
        d2 = H(j)*complex(d/sqrt(2),d/sqrt(2));
        d3 = H(j)*complex(0,d);
        d4 = H(j)*complex(-d/sqrt(2),d/sqrt(2));
        d5 = H(j)*complex(-d,0);
        d6 = H(j)*complex(-d/sqrt(2),-d/sqrt(2));
        d7 = H(j)*complex(0,-d);
        d8 = H(j)*complex(d/sqrt(2),-d/sqrt(2));
        
        p1 = abs(DFT_Y(j,i)-d1);
        p2 = abs(DFT_Y(j,i)-d2);
        p3 = abs(DFT_Y(j,i)-d3);
        p4 = abs(DFT_Y(j,i)-d4);
        p5 = abs(DFT_Y(j,i)-d5);
        p6 = abs(DFT_Y(j,i)-d6);        
        p7 = abs(DFT_Y(j,i)-d7);
        p8 = abs(DFT_Y(j,i)-d8);
        
        P = [p1 p2 p3 p4 p5 p6 p7 p8];
        [m,index] = min(P); 
        
        if (index==1)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 0;
        elseif (index==2)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 1;
        elseif (index==3)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 1; 
        elseif (index==4)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 0;
        elseif (index==5)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 0;
        elseif (index==6)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 1;
        elseif (index==7)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 1;
        elseif (index==8)
            rx_bitstream((i-1)*Nc*3+3*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*3+3*(j-1)+3) = 0;
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
ber = error/total;


end  
        
        

