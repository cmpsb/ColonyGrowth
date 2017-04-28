function [ipksi,AR] = corF(data,N)

fr = data(find(data(:,2)==N,1),1);
frames = find(data(:,1)==fr);

xpos = data(frames,3);
ypos = data(frames,4);
angle = data(frames,5);
AR = mean(data(frames,6));

D = ceil(max([max(xpos)-min(xpos) max(ypos)-min(ypos)])); %diameter of colony at framenr
stepsize = 1/5;
range = linspace(0,D,(D/stepsize));
ksi = zeros(length(xpos),length(range)-1);
ksi(:,1) = 1;

for i = 1:N
    c_x = xpos(i);
    c_y = ypos(i);
    c_angle = angle(i);
    xrelpos = xpos - c_x; %center colony around cell i
    yrelpos = ypos - c_y;
    dtc = sqrt(xrelpos.^2 + yrelpos.^2); % distance to cell
    for step = 2:(length(range)-1)
        oop = [];
        for nb = find(dtc > range(step) & dtc < (range(step+1)))
            dtheta = abs(c_angle - angle(nb));
            if dtheta > pi/2
                dtheta = pi - dtheta;
            end
            oop = [oop cos(2*dtheta)];
        end
        if ~isempty(oop)
            ksi(i,step) = mean(oop);
        else
            ksi(i,step) = ksi(i,step-1);
        end
    end
end

ksiMean = mean(ksi);
range = (range(1:end-1)+(stepsize/2))./D;
plot(range,ksiMean,':k')

ipx = linspace(0,0.2,100);
ipksi = interp1(range,ksiMean,ipx);
ipksi(isnan(ipksi))=1;

end


