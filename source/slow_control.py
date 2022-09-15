import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import datetime
import matplotlib.dates as mdates

class SlowControl():
    def __init__(self, path="/scratch/gpfs/aj9512/jadwin-365/data/temp/", selection="*"):
        self.path = path
        self.selection = selection
        self.files = glob.glob(self.path+self.selection)
    
    def read_temp_files(self):
        self.temp = []
        for i,f in enumerate(self.files):
            print(i,f)

            #read excel file into pandas dataframe
            df = pd.read_excel("{}".format(f), header=None, usecols=[0,1,2],skiprows=4)

            #get start time from each separate file
            start = pd.read_excel("{}".format(f), header=None, usecols=[1])
            start = start[1][1]
            start = datetime.datetime.strptime('{}'.format(start), '%a %b %d %H:%M:%S EDT %Y')

            #get time stamps in absolute time
            df[0] = [start + datetime.timedelta(minutes=x) for x in np.array(df[0]/1000.0/60.0)]

            #append to list of dataframes
            self.temp.append(df)

        #combine all dataframes into one
        self.temp = pd.concat(self.temp)

        # sort along timestamp axis
        self.temp = self.temp.sort_values(0,axis=0)
    
    def plot_temp(self,style='p'):
        fig, ax = plt.subplots()
        fig.autofmt_xdate()
        date_form = mdates.DateFormatter("%H:%M")
        ax.xaxis.set_major_formatter(date_form)

        if style=='p':
            plt.plot(self.temp[0], self.temp[1], label='Sensor 1')
            plt.plot(self.temp[0], self.temp[2], label='Sensor 2')
        elif style=='s':
            plt.scatter(self.temp[0], self.temp[1], label='Sensor 1', s=1)
            plt.scatter(self.temp[0], self.temp[2], label='Sensor 2', s=1)

        plt.hlines(xmin=np.min(self.temp[0]), xmax=np.max(self.temp[0]), y=87.3, linestyles='--', lw=1, color='k')

        plt.xlabel('Time [hh:mm]')
        plt.ylabel('Temperature [K]')
        plt.legend(loc='upper right')
        plt.show()
