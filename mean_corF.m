clear corm cor
figure;
hold on

data = {data0 data1 data2 data3 data4};
corm = cell(1,3); %store in cell for different N
loop = 0;
xplot = linspace(0,0.2,100); %change according to ipx in function corF

for nop = 2000
    loop = loop + 1;
    
    cor = zeros(length(data),length(xplot));
    for i = 1:length(data)
        [cor(i,:),AR] = corF(data{i},nop);
    end
    corm{loop} = mean(cor);
    error = std(cor);

    % p = polyfit(x,sopm(11,:),1);
    % fit = polyval(p,x);
    % plot(x,fit,'--b');

    plot(xplot,corm{loop},'-.','LineWidth',1.5);
    % x = 1:100:1000;
    % errorbar(x,nod(11,x),error(x));
    % axis([-100 max(xaxis)+100 -0.1 1.05]);
end

xlabel('Normalized distance to cell');
ylabel('Correlation in orientation');
grid on
hold off