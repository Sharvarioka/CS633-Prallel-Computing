import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

input_file="plot"+sys.argv[1]+".csv"

plot_df=pd.read_csv(input_file, error_bad_lines=False)

fig,ax=plt.subplots()
sns.set_style('ticks')
fig.set_size_inches(11.7, 8.27)
someNum = [i for i in range(0,105)]

color1 = ["#96cac1"]
color2= ["#5f9e6e"]
color3=["#cc8963"]


def plotTheBox(res, width, lineColor, lineWidth,colpal):
	sns.boxplot(y='Time',x='N',data=res,hue='Type',width=width,palette=sns.set_palette(sns.color_palette(colpal,1)))
	groupedResult=res.groupby('N')['Time'].median()
	ax.plot(groupedResult.values,color=lineColor,linewidth=lineWidth)


def getIndexToDrop(N):
	return list(filter(lambda x: x%3 != N , someNum))

# sns.boxplot(y='Time',x='N',data=plot_df,hue='Type')

res0=plot_df.drop(plot_df.index[getIndexToDrop(0)])
plotTheBox(res0, 0.4, '#96cac1', 1,color1)
# 'r-o'

res1=plot_df.drop(plot_df.index[getIndexToDrop(1)])
plotTheBox(res1, 0.2, '#5f9e6e', 1,color2)

res2=plot_df.drop(plot_df.index[getIndexToDrop(2)])
plotTheBox(res2, 0.15, "#cc8963", 1, color3)

plt.savefig("plot"+sys.argv[1]+".jpg")
plt.title('Box plot for '+sys.argv[1]+' processes')
plt.xlabel('Data points per process') 
plt.ylabel('log2 time in seconds') 
plt.savefig("plot"+sys.argv[1]+".jpg")
#plt.show()




