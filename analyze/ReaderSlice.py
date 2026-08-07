import glob
from Reader3D import Reader3D
import math

class SliceReader(Reader3D):

    def __init__(self, path, var_name: str, x_var: str, x_grid, y_var: str, y_grid, z_var: str, z_grid, color=None):
        self._files = glob.glob(path+'/'+var_name+'_slice_*')

        super(SliceReader, self).__init__(var_name, x_var, x_grid, y_var, y_grid, z_var, z_grid, color)


    def _grid_resolution(self):
        """Reads the name of the first file found to find information about grid resolution
        Returns [n_r, n_theta, n_time]
        """
        label = self._files[0].split('_slice_')[1]
        values = label.split("_")
        n_r = int(values[0])
        n_theta = int(values[1])
        n_time = int(len(self._files))

        print(f'Resolution {n_r, n_theta, n_time}')
        # t_digits = int(len(values[2]))
        # t_jump = 0

        # if len(files) > 1:
        #     label2 = files[1].split('_slice_')[1]
        #     values2 = label2.split("_")
        #     t_jump = int(values2[2]) - int(values[2])

        return [n_r, n_theta, n_time]
    

    def _find_time(self, lines):
        t = -1.
        for line in lines:
            if 'coord time' in line:
                t = float(line.split('coord time')[1].split()[0])
                break 
        return t


    def _find_value(self, lines, r_step, th_step):
        """Assuming grids index from 0
        Returns [r, theta, func_val]
        """
        for _ in range(100):
            if lines[0][0] == '#':
                lines.pop(0)
            else:
                break
            
        row = r_step * (self.resolution[1] + 1) + th_step
        th = (2 * th_step + 1)/(2 * self.resolution[1]) * math.pi
        vals = [float(x) for x in lines[row].split()]
        return [vals[0], th, vals[2]]