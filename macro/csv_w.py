import csv

with open('/Users/garidor/Desktop/first-year/macro/capital_stock.csv', 'r') as f:
	reader = csv.reader(f, delimiter=',')
	with open('/Users/garidor/Desktop/first-year/macro/capital_stock_2.csv', 'wb') as wr:
		file_writer = csv.writer(wr)
		for row in reader:
			try:
				file_writer.writerow([float(row[1])])
				file_writer.writerow([float(row[1])])
				file_writer.writerow([float(row[1])])
				file_writer.writerow([float(row[1])])
			except:
				continue