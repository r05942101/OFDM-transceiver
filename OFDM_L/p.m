
bpsk = zeros(1,40);


for i =1:40
    bpsk(i) = OFDM_BPSK_L_Changed(i);
end

figure(1)
x = 1:40
plot(x,bpsk,'--x');
xlabel('CP length');
ylabel('BER');
grid on


