import numpy as np
from matplotlib import pyplot as plt



PART = 0
ITER = 0
xlims = []
ylims = []

with open("./_runresult.txt", 'r') as f:
    line = f.readline()
    line = f.readline()
    xlims.append(float(line.split()[0]))
    xlims.append(float(line.split()[1]))
    line = f.readline()
    ylims.append(float(line.split()[0]))
    ylims.append(float(line.split()[1]))
    line = f.readline()
    PART = int(line.split()[3])
    line = f.readline()
    ITER = int(line.split()[3])

    
with open("./_dump.txt", 'r') as f:
    def animate(i):
        xs = []
        ys = []
        for _ in range(PART):
            linee = f.readline()
            xs.append(float(linee.split()[0]))
            ys.append(float(linee.split()[1]))

        fig = plt.figure(figsize=(6, 6))
        plt.plot(xs, ys, 'b.' )
        plt.grid()
        plt.xlim(1.1 * xlims[0], 1.1 * xlims[1])
        plt.ylim(1.1 * ylims[0], 1.1 * ylims[1])
        plt.xlabel('x', fontsize = 14)
        plt.ylabel('y', fontsize = 14)
        plt.savefig(f'./img/img_{i}.png',
                    transparent = False,  
                    facecolor = 'white'
                )
        plt.close()
        return plt.savefig

    for i in range(ITER):
        animate(i)
