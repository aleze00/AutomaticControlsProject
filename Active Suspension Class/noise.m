function nu = noise(t)

g = 9.81; % [m/s^2] grevitational acceleration
std_pot = 0.001; % [m] potentiomenter standard deviation
std_acc = 0.05*g; % [m/s^2] accelerometer standard deviation

nu = [std_pot*randn(1); 
      std_acc*randn(1)];
end