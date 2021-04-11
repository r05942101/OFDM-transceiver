%{
bpsk = zeros(1,31);
parr_bpsk = zeros(1,31);
bapsk = zeros(1,31);
parr_8psk = zeros(1,31);
qam = zeros(1,31);
parr_qam = zeros(1,31);
for i =0:30
    parr_bpsk(i+1) = OFDM_BPSK2(i);
    parr_8psk(i+1) = ofdm_8psk2(i);
    parr_qam(i+1) = ofdm_16QAM2(i);
end
 %}
figure(1)
x = 0:30
plot(x,parr_bpsk,'--x',x,parr_8psk,'--x',x,parr_qam,'--x');
xlabel('E_{b}/N_{0}');
ylabel('PAPR');
legend('bpsk','8psk','16QAM');
grid on

