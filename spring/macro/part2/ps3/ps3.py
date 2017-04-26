import pandas
import heapq
import matplotlib.pyplot as plt

def mean(L):
	return(sum(L) / float(len(L)))

def calc_equity_premium(p, q, sigma, rho, theta, b):
	first = sigma*sigma * theta
	expect = mean(map(lambda x: (1 - x)**(-theta), b))
	second = p * (1 - q) * expect
	third = mean(map(lambda x: (1 - x)**(1-theta), b))
	fourth = mean(b)
	return first + second - third - fourth


df = pandas.read_excel('/Users/garidor/Desktop/first-year/spring/macro/part2/ps3/DataForPS3.xlsx')
col = "% fall in real per capita GDP"
b = df[col].values
p = 0.017
q = 0.4
gamma = 0.025
sigma = 0.02
rho = 0.03
theta = 4
vals = []
print(calc_equity_premium(p, q, sigma, rho, theta, b))
for theta in range(1, 7):
	vals.append(calc_equity_premium(p, q, sigma, rho, theta, b))
plt.plot([1, 2, 3, 4, 5, 6], vals)
#plt.show()

b = heapq.nsmallest((len(b) - 6), b)
print(calc_equity_premium(p, q, sigma, rho, theta, b))
b = heapq.nlargest((len(b) - 30), b)
print(calc_equity_premium(p, q, sigma, rho, theta, b))