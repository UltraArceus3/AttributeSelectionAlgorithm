
from matplotlib import pyplot as plt





def generatePRocessedData():
    # x-axis values
    x = [28637,35816,42515,53013]
    
    # Y-axis values
    y = [5.91324,8.94686,12.6739,19.5101]
    
    # Function to plot
    plt.plot(x,y)
    plt.ylabel('Time')
    plt.xlabel('# Records')
    
    # function to show the plot
    plt.show()


