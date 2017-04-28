function [sop,xn] = find_sopnew(datafile)

clear xn oop

timestep = 0;
timeTheta = 0;
aveTheta = 0;
sumTheta = 0;
N = 0;

for i = 1:(length(datafile(:,1))-1)
    timeTheta = timeTheta + datafile(i,5);
    dTheta = abs(aveTheta - datafile(i,5));
    if dTheta > (pi/2)
        dTheta = pi - dTheta;
    end
    sumTheta = sumTheta + (2 * cos(dTheta)^2 - 1);
    N = N + 1;
    if datafile(i,1) ~= datafile(i+1,1)
        timestep = timestep + 1;
        aveTheta = timeTheta / N;
        sop(timestep) = sumTheta / N;
        timeTheta = 0;
        sumTheta = 0;
        xn(timestep) = N;
        N = 0;
    end
end

%plot(xn,oop);

end