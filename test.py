import os
p = '~'
os.makedirs(p +"/omicDC_results", exist_ok=True)
file_name = 'my_file.txt'
f = open(p +"/omicDC_results/"+file_name, 'a+')  # open file in append mode
f.write('python rules')
f.close()