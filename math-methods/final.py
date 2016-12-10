import scipy.sparse as sps
import scipy.sparse.linalg as lg
import math

N = 1001
rows, cols = N, N
sps_acc = sps.lil_matrix((rows, cols))
for i in xrange(N):
	if i <= 10:
		if i == 0:
			sps_acc[0, 0] = 0.67
		else:
			sps_acc[i, i] = 0.33
			sps_acc[i, 0] = 0.34
		sps_acc[i, i+10] = 0.33
	elif i > 10 and i < 990:
		sps_acc[i, 0] = 0.01
		sps_acc[i, i - 10] = 0.33
		sps_acc[i, i] = 0.33
		sps_acc[i, i + 10] = 0.33
	elif i >= 990:
		sps_acc[i, 0] = 0.01
		sps_acc[i, i - 10] = 0.33
		if i == N-1:
			sps_acc[N-1, N-1] = float(0.66)
		else:
			sps_acc[i, i] = 0.33
			sps_acc[i, N-1] = 0.33

distriubtion = sps_acc**10000000

sum_of_vec = 0.0
var_of_vec = 0.0
skewness_of_vec = 0.0
for i in xrange(N):
	sum_of_vec += distribution[1, i]

mean = sum_of_vec / N
for i in xrange(N):
	var_of_vec += (distribution[1, i] - mean)**2

var_of_vec = var_of_vec / (N-1)
#std = math.sqrt(var_of_vec)
#for i in xrange(N):
#	skewness_of_vec += ((distribution[1, i] - mean)/std)^3


print(distribution)
print(mean)
print(var_of_vec)
