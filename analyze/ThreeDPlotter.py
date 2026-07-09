import matplotlib.pyplot as plt

class ThreeDPlotter:

    def __init__(self):
        self.var_name = ''
        self.x_var = ''
        self.y_var = ''
        self.y_grid = []
        self.z_var = ''

        self.xlimits = None
        self.ylimits = None
        self.legend_prefix = ''

        self.x_val = []
        self.y_val = []
        self.z_val = []

        self.num_plot = [[[]]]
        self.an_plot = [[[]]]
        self.aux_plot = [[[]]]

        self._delay = 0.01

    def set_xlimits(self, limits):
        self.xlimits = limits

    def set_ylimits(self, limits):
        self.ylimits = limits
    
    def set_legend_prefix(self, prefix):
        self.legend_prefix = prefix

    def set_delay(self, delay):
        self._delay = delay

    def plot_z(self, z, ax, numerical, analytical, auxiliary):
        """Plots information at a set z value given an pyplot axis
        Does not actually display the plot
        """

        for y_iteration in range(len(self.y_grid)):
            if len(self.y_grid) > 1:
                label = f'{self.legend_prefix}{self.y_var}={round(self.y_grid[y_iteration], 2)}'
            else:
                label = self.legend_prefix

            if numerical:
                ax.scatter(self.x_val, self.num_plot[z][y_iteration], label=f'num {label}', marker='.')

            if analytical:
                ax.plot(self.x_val, self.an_plot[z][y_iteration], label=f'an {label}')
                    
            if auxiliary:
                ax.scatter(self.x_val, self.aux_plot[z][y_iteration], label=f'aux {label}', marker='x')

        ax.set(xlabel=f'{self.x_var}', ylabel=f'{self.var_name}')
        ax.set_xlim(self.xlimits)
        ax.set_ylim(self.ylimits)
        plt.legend()
        plt.title(f"{self.z_var} = {round(self.z_val[z], 2)}")
    
    def plot(self, numerical=True, analytical=False, auxiliary=False):
        "Runs animated plots based off the values found"
        _, ax = plt.subplots(1)
        for z_iteration in range(len(self.z_grid)): #type: ignore
            ax.clear()
            self.plot_z(z_iteration, ax, numerical, analytical, auxiliary)
            plt.pause(self._delay)

        plt.show()