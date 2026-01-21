function IncrementalScatter(pts,turns)
hold on;
n = size(pts,1);
for i=0:turns;
    scatter(ordered(1+i*n/turns:(i+1)*n/turns,1), ordered(1+i*n/turns:(i+1)*n/turns,2));
    pause;
end
hold off;
end