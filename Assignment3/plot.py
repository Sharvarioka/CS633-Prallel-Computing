import matplotlib.pyplot as plt 
    
# x axis values 
x = ['1_1','1_2','1_4','2_1','2_2','2_4'] 

f = open('output.txt')

y_coord = []

read_lines = f.readlines()
i=1
# corresponding y axis values 
for line in read_lines:
	if i%3==0:
		y_coord.append(float(line.strip()))
	i=i+1

f.close()
plt.plot(x, y_coord) 
plt.xlabel('Node_cores per node') 
plt.ylabel('Time in seconds') 
plt.title('Plot for time comparison') 
plt.savefig('output.jpg')
plt.show() 