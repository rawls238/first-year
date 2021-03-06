import math
import copy
from decimal import *
import matplotlib.pyplot as plt

def x_n(a, b, x_n_1, y_n_1):
	return a - x_n_1**2 + b * y_n_1

def y_n(x_n_1):
	return x_n_1

def calculate_dynamics(x_0, y_0, num_iter):
	current_x = Decimal(x_0)
	current_y = Decimal(y_0)
	a = Decimal(1.4)
	b = Decimal(0.3)
	x_vals = [current_x]
	y_vals = [current_y]
	for i in range(1, num_iter):
		prev_y = current_y
		prev_x = copy.copy(current_x)
		current_x = x_n(a, b, prev_x, prev_y)
		current_y = x_vals[i - 1]
		x_vals += [current_x]
		y_vals += [current_y]
	#line, = plt.plot(range(0, num_iter+1), x_vals)
	#plt.show()
	#line, = plt.plot(x_vals, y_vals)
	#plt.show()
	#plt.hist(x_vals, bins=40)
	#plt.show()
	return x_vals, y_vals


getcontext().prec = 57
x12_val, y12_val = calculate_dynamics(0.0, 0.0, 5000)
x2_val, y2_val = calculate_dynamics(10.0**-8, 10.0**-8, 5000)
diff = [math.log(abs(x2_val[i] - x1_val[i])) for i in range(0, len(x2_val))]
plt.plot(range(0,5001), diff)
plt.show()