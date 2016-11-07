import csv

with open('/Users/garidor/Desktop/randomwage.csv', 'r') as f:
	reader = csv.reader(f, delimiter=',')
	with open('/Users/garidor/Desktop/randomwage2.csv', 'wb') as wr:
		file_writer = csv.writer(wr)
		for row in reader:
			for item in row:
				file_writer.writerow([item])