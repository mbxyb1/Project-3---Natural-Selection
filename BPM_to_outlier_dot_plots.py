# BPM_to_outlier_dot_plots
# Move everything you need into a 'BPM_to_outlier_dot_plots' folder and run inside it:
	# It will contain a 'temp_files' folder
	# 
# Move your input BPM from ScanTools into 'BPM_to_outlier_dot_plots' and rename it 'input_BPM.txt'

import argparse, os, subprocess

parser=argparse.ArgumentParser(description="Extracts the top percentage of outliers based on one or more metrics.")
parser.add_argument('-n', type=str, metavar='run_name', required=True, help='A name for the run e.g. Rho or Fst')
#parser.add_argument('-c', type=str, nargs='+', metavar='metric_columns', required=True, help='The column(s) which contain the metrics that you want the top outliers from. Location is zero-based i.e. the first column is "0", the second "1", the third "2", etc')
#parser.add_argument('-t', type=str, nargs='+', metavar='tops_bottoms', required=True, help='For each column, in the order you entered the columns, write "Top" for highest outliers or "Bottom" for lowest outliers.')
parser.add_argument('-p', type=str, metavar='percentage', required=True, help='The percentage (0-100) of outliers you want as a fraction, i.e. 1 will give you the genes which are in the top 1 percent outliers in all the metrics specified with -c.')
#parser.add_argument('-x', type=str, metavar='exclusive', required=False, default='True', help='For an exclusive list set as True (only sites that are outliers for all metrics), for an inclusive list set to False (sites that are an outlier in one or more metrics).')
args=parser.parse_args()

#c = ' '.join(args.c)
#t = ' '.join(args.t)

#os.system('python3 Find_BPM_outliers.py -i input_BPM.txt -c '+c+' -t '+t+' -p '+args.p+' -x '+args.x+' -o temp_files/temp1.txt') #list of outlier windows
os.system('bedtools intersect -a temp_files/LABvsNENtemp.txt -b cp_genes_only.gft -wb | cut -d = -f 2 > temp_files/temp2.txt') #list of outlier genes #modified already 

#remove duplicates
file_to_sort=open('temp_files/temp2.txt', 'r')
list_to_check=[]
for count4, gene2 in enumerate(file_to_sort):
	list_to_check.append(gene2)
file_to_sort.close()
reduced_list=list(dict.fromkeys(list_to_check))
temp3=open('temp_files/temp3.txt', 'w+')
for count5, gene3 in enumerate(reduced_list):
	temp3.write(gene3)
temp3.close()

#Sort gene list
temp3=open('temp_files/temp3.txt', 'r')
temp4=open('temp_files/temp4.txt', 'w+')
gene_list=[]
for count2, line in enumerate(temp3):
	gene_list.append(int(line.replace('g', '')))
gene_list.sort()
for count3, item in enumerate(gene_list):
	temp4.write('g'+str(item)+'\n')
temp3.close()
temp4.close()

os.system('mkdir 0_'+args.n+'_'+args.p+'percent') #####creates that directory 

#grep lines for exact matches one by one (to return all lines for lists that contain e.g. g1 and g10, g10 is missed with grep -w)

temp4=open('temp_files/temp4.txt', 'r')
### Annotated_out=open('0_'+args.n+'_'+args.p+'percent/Annotated_genes_'+args.n+'_'+args.p, 'w+')
temp5=open('temp_files/temp5.txt', 'w+')

for count6, gene4 in enumerate(temp4):
	gene4=gene4.replace('\n', '')
	### result = subprocess.run(['grep', '-w', '%s'%gene4, 'Annotated_Hits.txt'], stdout=subprocess.PIPE)
	### Annotated_out.write(result.stdout.decode())
	result = subprocess.run(['grep', '-w', '%s'%gene4, 'cp_genes_only.gft'], stdout=subprocess.PIPE)
	temp5.write(str(result.stdout.decode()))
	
temp4.close()
### Annotated_out.close()
temp5.close()

#os.system('grep -w -f temp_files/temp4.txt Annotated_Hits.txt > 0_'+args.n+'_'+args.p+'percent/Annotated_genes_'+args.n+'_'+args.p)
#os.system('grep -w -f temp_files/temp4.txt genes_only.gtf > temp_files/temp5.txt')

os.system('cut -f 1,4-5 temp_files/temp5.txt > temp_files/temp6.txt') #list of gene locations

os.system('paste temp_files/temp4.txt temp_files/temp4.txt temp_files/temp6.txt > candidate_gene_file.txt')

#Change this bit as needed
#os.system('python3 R.gene_window_plot.py tetraploids 96 diploids 16 gene_orientation_file.txt candidate_gene_file.txt 30000') ###this worked fine 

os.system('python3 R.gene_window_plot2.py LAB 16 NEN 1 cp_gene_orientation_file.txt candidate_gene_file.txt 100')

os.system('mv *.pdf 0_'+args.n+'_'+args.p+'percent/')
os.system('mv *.out 0_'+args.n+'_'+args.p+'percent/')
os.system('mv candidate_gene_file.txt 0_'+args.n+'_'+args.p+'percent/')
#os.system('rm *.R')



