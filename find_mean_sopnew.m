clear xaxis L sopmean

[sop,nop] = find_sopnew(data);
[sop1,nop1] = find_sopnew(data1);
[sop2,nop2] = find_sopnew(data2);
[sop3,nop3] = find_sopnew(data3);
[sop4,nop4] = find_sopnew(data4);

L = min([length(sop),length(sop1),length(sop2),length(sop3),length(sop4)]);
if length(sop) == L
    xaxis = nop;
elseif length(sop1) == L
    xaxis = nop1;
elseif length(sop2) == L
    xaxis = nop2;
elseif length(sop3) == L
    xaxis = nop3;
else
    xaxis = nop4;
end

sopmean = zeros(6,L);
sopmean(1,:) = sop(1:L);
sopmean(2,:) = sop1(1:L);
sopmean(3,:) = sop2(1:L);
sopmean(4,:) = sop3(1:L);
sopmean(5,:) = sop4(1:L);

for i = 1:L
    sopmean(6,i) = mean(sopmean(1:5,i));
end

figure;
plot(xaxis,sopmean(6,:),'x');
axis([-50 max(xaxis)+100 -0.1 1.05]);
ylabel('Scalar Order Parameter');
xlabel('Number of particles');