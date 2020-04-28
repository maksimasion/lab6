import matplotlib.pyplot as plt
import pylab
import numpy as np

def makeData():
    X = np.arange(-40, 25, 0.1)
    Y = np.arange(-30, 15, 0.1)
    X, Y = np.meshgrid(X, Y)

    Z = X**4 + 20.*(X**3) + 2*(X**2)*(Y**2) + 36*(X**2)*Y + 312*(X**2)+ 20*X*(Y**2) + 360*X*Y + 2121*X + Y**4+ 36*(Y**3)+ 537*(Y**2) + 3834*Y + 11308
    
    return X, Y, Z


if __name__ == '__main__':
    x, y, z = makeData()
    fig = pylab.figure()

    data=np.loadtxt ("tr.txt")
    X=data[:,0]
    Y=data[:,1]
    k=len(X)
    X1=[X[k-1]]
    Y1=[Y[k-1]]

    plt.plot(X,Y,"g", linewidth = 2)
    plt.plot(X1,Y1,"og")
    plt.xlabel('ось x1')
    plt.ylabel('ось x2')

    selected_x = X1[0]
    selected_y = Y1[0]

    pylab.rc('font', family = 'verdana')

    arrowprops = {
        'arrowstyle': '->',
    }

    a = np.array ([X[k-1], Y[k-1]-6])

    pylab.annotate (u'Точка минимума',
                    xy=(X1[0], Y1[0]),
                    xytext = a,
                    arrowprops = arrowprops)
    fig.set_figwidth(8)
    fig.set_figheight(5)
    cs = pylab.contour(x, y, z, 50)
    pylab.clabel(cs)
    pylab.show()
