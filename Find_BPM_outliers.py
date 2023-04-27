#Extracts the top outliers from a BPM list generated by ScanTools

import argparse

parser=argparse.ArgumentParser(description="Extracts the top percentage of outliers based on one or more metrics.")
parser.add_argument('-i', type=str, metavar='input_BPM', required=True, help='Input between population metrics (BPM) file from ScanTools. File is tab deliminated (e.g. outname	scaff	start	end	win_size	num_sites	num_snps	Rho	FstWC	dxy	AFD	FixedDiff	FstH	FstN.')
parser.add_argument('-c', type=int, nargs='+', metavar='metric_columns', required=True, help='The column(s) which contain the metrics that you want the top outliers from. Location is zero-based i.e. the first column is "0", the second "1", the third "2", etc')
parser.add_argument('-t', type=str, nargs='+', metavar='tops_bottoms', required=True, help='For each column, in the order you entered the columns, write "Top" for highest outliers or "Bottom" for lowest outliers.')
parser.add_argument('-p', type=float, metavar='percentage', required=True, help='The percentage (0-100) of outliers you want as a fraction, i.e. 1 will give you the genes which are in the top 1 percent outliers in all the metrics specified with -c.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='A name and location for the output file.')
parser.add_argument('-x', type=bool, metavar='exclusive', required=False, default=True, help='For an exclusive list set as True (only sites that are outliers for all metrics), for an inclusive list set to False (sites that are an outlier in one or more metrics).')
args=parser.parse_args()

output_file=open(args.o, 'w+')

columns=list(args.c)
top_or_bottom=list(args.t)
all_locations=[]
all_top_outliers=[]
to_print=[]

for count, column in enumerate(columns):
	input_file=open(args.i, 'r')
	values=[]
	locations=[]
	top_outliers=[]
	for count1, line in enumerate(input_file):
		if count1 != 0: #Escape the headder
			line=line.replace('\n', '')
			line=line.split('\t')
			values.append(line[column])
			locations.append(line[1]+'\t'+line[2]+'\t'+line[3]+'\n')
	sorted_values=list(dict.fromkeys(values))
	if top_or_bottom[count] == 'Top':
		sorted_values.sort(reverse=True)
	if top_or_bottom[count] == 'Bottom':
		sorted_values.sort()
	if top_or_bottom[count] != 'Top' and top_or_bottom[count] != 'Bottom':
		print('-t not specified correctly, you must put "Top" or "Bottom" seperated by a space and without the quotes, -t is case sensative')
	number_of_outliers=round((len(sorted_values)/100)*args.p)
	for count2, value in enumerate(sorted_values):
		if count2 <= number_of_outliers:
			top_outliers.append(value)
	for count3, metric in enumerate(top_outliers):
		while metric in values:
			index=values.index(metric)
			all_top_outliers.append(locations[index])
			del locations[index]
			del values[index]
	input_file.close()

all_top_outliers.sort()

if args.x == True:
	for count4, location in enumerate(all_top_outliers):
		if all_top_outliers.count(location) == len(columns) and location not in to_print:
			to_print.append(location)

if args.x == False:
	to_print=list(dict.fromkeys(all_top_outliers))

for count5, thing in enumerate(to_print):
	output_file.write(thing.replace('.0', ''))

output_file.close()