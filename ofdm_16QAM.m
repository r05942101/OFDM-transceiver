function ber = ofdm_16QAM(Eb_N0)
load LTI
% define parameter
N0 = 2;
Nc = 128;
L = 30;
num = 100;
total = num*Nc*4;

s_bitstream = randi([0,1],1,total);

p_bistream = zeros(Nc*4,num); 
% serial to parallel
for i =1:num
    for j =1:Nc*4
        p_bistream(j,i) = s_bitstream(Nc*4*(i-1)+j);
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
row = size(p_bistream);
for i =1:num
% bpsk modulation
    U = zeros(1,Nc);
    for j =1:4:(row-3)
        if(p_bistream(j,i)==0 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(-3*d,-3*d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(-3*d,-d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(-3*d,d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(-3*d,3*d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(-d,-3*d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(-d,-d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(-d,d);
        elseif (p_bistream(j,i)==0 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(-d,3*d);
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(d,-3*d);   
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(d,-d);        
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(d,d);        
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==1 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(d,3*d);        
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(3*d,-3*d);   
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==0  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(3*d,-d); 
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==1)
            U((j-1)/4+1)=complex(3*d,d);      
        elseif (p_bistream(j,i)==1 && p_bistream(j+1,i)==0 && p_bistream(j+2,i)==1  && p_bistream(j+3,i)==0)
            U((j-1)/4+1)=complex(3*d,3*d);
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
        a1 = H(j)*complex(-3*d,-3*d);
        a2 = H(j)*complex(-3*d,-d);
        a3 = H(j)*complex(-3*d,d);
        a4 = H(j)*complex(-3*d,3*d);
        a5 = H(j)*complex(-d,-3*d);
        a6 = H(j)*complex(-d,-d);
        a7 = H(j)*complex(-d,d);
        a8 = H(j)*complex(-d,3*d);
        a9 = H(j)*complex(d,-3*d);
        a10 = H(j)*complex(d,-d);
        a11 = H(j)*complex(d,d);
        a12 = H(j)*complex(d,3*d);
        a13 = H(j)*complex(3*d,-3*d);
        a14 = H(j)*complex(3*d,-d);
        a15 = H(j)*complex(3*d,d);
        a16 = H(j)*complex(3*d,3*d);

        d1=abs(DFT_Y(j,i)-a1);
        d2=abs(DFT_Y(j,i)-a2);
        d3=abs(DFT_Y(j,i)-a3);
        d4=abs(DFT_Y(j,i)-a4);
        d5=abs(DFT_Y(j,i)-a5);
        d6=abs(DFT_Y(j,i)-a6);        
        d7=abs(DFT_Y(j,i)-a7);
        d8=abs(DFT_Y(j,i)-a8);
        d9=abs(DFT_Y(j,i)-a9);
        d10=abs(DFT_Y(j,i)-a10);
        d11=abs(DFT_Y(j,i)-a11);
        d12=abs(DFT_Y(j,i)-a12);
        d13=abs(DFT_Y(j,i)-a13);
        d14=abs(DFT_Y(j,i)-a14);
        d15=abs(DFT_Y(j,i)-a15);
        d16=abs(DFT_Y(j,i)-a16);
        

        D = [d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16];
        [m,index] = min(D);
     
        if (index==1)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
        elseif (index==2)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        elseif (index==3)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1; 
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        elseif (index==4)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
        elseif (index==5)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
        elseif (index==6)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        elseif (index==7)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        elseif (index==8)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
        elseif (index==9)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
        elseif (index==10)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        elseif (index==11)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
         elseif (index==12)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
         elseif (index==13)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
         elseif (index==14)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
         elseif (index==15)
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 1;
        else 
            rx_bitstream((i-1)*Nc*4+4*(j-1)+1) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+2) = 0;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+3) = 1;
            rx_bitstream((i-1)*Nc*4+4*(j-1)+4) = 0;
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
        
        

