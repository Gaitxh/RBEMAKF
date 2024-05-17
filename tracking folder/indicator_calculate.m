function [pos, vel] = indicator_calculate(sim,res)
pos = (sim.x(1,:)-res.x(1,:)).^2+(sim.x(2,:)-res.x(2,:)).^2;
vel = (sim.x(3,:)-res.x(3,:)).^2+(sim.x(4,:)-res.x(4,:)).^2;
end