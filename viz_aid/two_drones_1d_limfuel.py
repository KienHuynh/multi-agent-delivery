import numpy as np
import matplotlib.pyplot as plt


def find_intersection(xy1, xy2):
    m1 = (xy1[0,1]-xy1[1,1])/(xy1[0,0]-xy1[1,0])
    b1 = xy1[0,1]-m1*xy1[0,0]
    
    m2 = (xy2[0,1]-xy2[1,1])/(xy2[0,0]-xy2[1,0])
    b2 = xy2[0,1]-m2*xy2[0,0]
    
    x = (b1 - b2) / (-m1 + m2)
    y = (b1*m2 - b2*m1) / (-m1 + m2)
    
    return np.asarray([x, y])
    

def fig1(xy1, xy2):
    # Draw the line segments
    plt.figure(1)
    plt.axis('equal')
    d1_plot = plt.plot(xy1[:,0], xy1[:,1], color='b')
    d2_plot = plt.plot(xy2[:,0], xy2[:,1], color=[1, 0.5, 0])
    
    # Draw the vertical axis
    plt.plot([0,0], [0,10], color='k')
    
    # Draw the starting line
    start_plot = plt.plot([0,8], [2,2], color='r')
    
    plt.legend((d1_plot[0], d2_plot[0], start_plot[0]), ('Drone 1', 'Drone 2', 'Package'))
    

def fig2(xy1, xy2, xy1_r):
    # Find the intersection after drone 1 goes back up
    xy_int = find_intersection(xy1_r, xy2)
    
    # Draw the line segments
    plt.figure(2)
    plt.axis('equal')
    d1_plot = plt.plot(xy1[:,0], xy1[:,1], color='b')
    d2_plot = plt.plot(xy2[:,0], xy2[:,1], color=[1,0.5,0])
    plt.plot(xy1_r[:,0], xy1_r[:,1], color='b')
    
    # Draw the vertical axis
    plt.plot([0,0], [0,10], color='k')
    
    # Draw the starting line
    start_plot = plt.plot([0,8], [2,2], color='r')
    
    # Draw the target lines
    end_plot = plt.plot([0,8], [xy1[1,1], xy1[1,1]], color='g')
    plt.plot([0,8], [xy_int[1], xy_int[1]], color='g')
    
    plt.legend((d1_plot[0], d2_plot[0], start_plot[0], end_plot[0]), ('Drone 1', 'Drone 2', 'Package', 'Target'))

    
def fig3(xy1, xy2, xy1_r):
    # Find the intersection after drone 1 goes back up
    xy_int = find_intersection(xy1_r, xy2)
    
    #xy1 = [[xy1[0,0], [xy_int[0,1]], []]
    
    # Draw the line segments
    plt.figure(3)
    plt.axis('equal')
    d1_plot = plt.plot(xy1[:,0], xy1[:,1], color='b')
    d2_plot = plt.plot(xy2[:,0], xy2[:,1], color=[1,0.5,0])
    plt.plot(xy1_r[:,0], xy1_r[:,1], color='b')
    
    # Draw the vertical axis
    plt.plot([0,0], [0,10], color='k')
    
    # Draw the starting line
    start_plot = plt.plot([0,8], [2,2], color='r')
    
    # Draw the target lines
    end_plot = plt.plot([0,8], [xy1[1,1], xy1[1,1]], color='g')
    plt.plot([0,8], [xy_int[1], xy_int[1]], color='g')
    
    plt.legend((d1_plot[0], d2_plot[0], start_plot[0], end_plot[0]), ('Drone 1', 'Drone 2', 'Package', 'Target'))
    

if __name__ == '__main__':
    # Max range 3
    # Speed 1
    xy1 = np.asarray([[0, 4], [4, 0]])
    # Drone 1 after it reflects on the starting line
    xy1_r = np.asarray([[2, 2], [4, 4]])
    
    # Max range 6
    # Speed 2.2
    xy2 = np.asarray([[0, 8], [3, 1.4]])
    
    plt.ion()
    fig1(xy1, xy2)
    fig2(xy1, xy2, xy1_r)
    fig3(xy1, xy2, xy1_r)
    plt.show()
    plt.pause(0)