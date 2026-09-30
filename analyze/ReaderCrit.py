import glob
import matplotlib.pyplot as plt


class CritReader():

    def __init__(self, path: str, file_extension: str, x_var: str, x_col: int, y_var: str, y_col: int):
        self._files = glob.glob(path+'/**/output/*'+file_extension)
        
        self.x_col = x_col
        self.y_col = y_col

        self.x_var = x_var.lower()
        self.legend_prefix = ''
        self.y_var = y_var.lower()

        self.xlimits = None
        self.ylimits = None

        self.x_val = [[]] # ie. The points on x axis
        self.y_val = [] # ie. The legend for each line
        
        self.num_plot = [[]]
        # self.an_plot = [[]]
        # self.aux_plot = [[]]
    
    def read(self):
        self.x_val = [[] for _ in range(len(self._files))] # type: ignore

        self.num_plot = [[] for _ in range(len(self._files))] # type: ignore
        # self.an_plot = [[] for _ in range(len(self._files))] # type: ignore
        # self.aux_plot = [[] for _ in range(len(self._files))] # type: ignore

        for i, file in enumerate(self._files):
            print(f"Reading {file}")
            f = open(file)
            lines = f.readlines()

            eta = file.split('/output/')[0].split('/')[-1]
            print(eta)
            self.y_val.append(eta)

            for line in lines:
                if line[0] == '#':
                    continue

                vals = [float(x) for x in line.split()]
                self.x_val[i].append(vals[self.x_col])

                self.num_plot[i].append(vals[self.y_col])
                
    def set_xlimits(self, limits):
        self.xlimits = limits


    def set_ylimits(self, limits):
        self.ylimits = limits

    
    def set_legend_prefix(self, prefix):
        self.legend_prefix = prefix

    
    def plot(self):
        "Plots all y vals onto the x value"
        _, ax = plt.subplots(1)
        for y in range(len(self.y_val)):
            if len(self.num_plot) > 1:
                label=f'{self.legend_prefix} {self.y_val[y]}'
            else:
                label=self.legend_prefix
                
            ax.plot(self.x_val[y], self.num_plot[y], label=f'{label}')
                
            # if analytical:
            #     ax.plot(self.x_val[y], self.an_plot[y], label=f'an {label}')
                
            # if auxiliary:
            #     ax.scatter(self.x_val[y], self.aux_plot[y], label=f'aux {label}')
        ax.set_xlim(self.xlimits)
        ax.set_ylim(self.ylimits)
