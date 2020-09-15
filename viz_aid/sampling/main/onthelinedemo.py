import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


def solve_quad(xp, yp, s, t, v0, v1):
    xq1 = (v0 * (v0 * xp + (- v0 ** 2 * yp ** 2 + v1 ** 2 * xp ** 2 + v1 ** 2 * yp ** 2) ** (1. / 2.))) / (
                v0 ** 2 - v1 ** 2)
    xq2 = (v0 * (v0 * xp - (- v0 ** 2 * yp ** 2 + v1 ** 2 * xp ** 2 + v1 ** 2 * yp ** 2) ** (1. / 2.))) / (
                v0 ** 2 - v1 ** 2)
    if t > xq1 > s:
        return xq1
    if t > xq2 > s:
        return xq2
    return -1


def nat_heu(v0, v1, delay):
    # handoff calculation
    R = np.sqrt(p1[0] ** 2 + p1[1] ** 2)
    phi = np.arccos(p1[0] / R)
    c = delay * v0
    theta = np.linspace(0, 2 * 3.14156, 500)
    cos_ = np.cos(phi - theta)

    r = (c*(v1**2) -
         v0*((R**2)*(v1**2) -
             (R**2)*(v0**2) +
             (c**2)*(v1**2) +
             (R**2)*(v0**2)*(cos_**2) +
             2*R*c*(v1**2)*cos_)**(1/2) +
         R*(v0**2)*cos_)/(v0**2 - v1**2)

    # a = (1 / (v1 ** 2) - 1 / (v0 ** 2))
    # b = -2 * (R * cos_ / (v1 ** 2) - delay/v0)
    # c = (R**2) / (v1 ** 2) - (delay ** 2)
    # delta = (b ** 2) - 4 * a * c;
    # #r1 = (-b + sqrt(delta)) / (2 * a);
    # r = (-b + np.sqrt(delta)) / (2 * a)

    qx = r * np.cos(theta)
    qy = r * np.sin(theta)
    makespans = (
                        np.sqrt((p1[0] - qx) ** 2 + (p1[1] - qy) ** 2) +
                        np.sqrt((t[0] - qx) ** 2 + qy ** 2)
                ) / v1
    a0_t = delay + t[0] / v0
    a1_t = (
                        np.sqrt((p1[0] - s[0]) ** 2 + (p1[1] - s[1]) ** 2) +
                        t[0]
                ) / v1

    best_ind = np.argmin(makespans)
    best_makespan = makespans[best_ind]
    if best_makespan > a0_t:
        if (a1_t > a0_t):
            return []
        else:
            return []
    q1 = np.asarray([qx[best_ind], qy[best_ind]])

    # Compute natural heuristic
    #xq1 = solve_quad(p1[0], p1[1], s[0], t[0], v0, v1)
    #if (xq1 == -1):
    #    return []

    nat_q1 = np.asarray([r[0], 0])
    nat_t = (
                    np.sqrt(np.sum((nat_q1 - p1) ** 2)) +
                    np.sqrt(np.sum((nat_q1 - t) ** 2))
            ) / v1
    approx_ratio = nat_t / best_makespan

    return [q1, nat_q1, best_makespan, nat_t, approx_ratio, qx, qy]

def update(val):
    v0 = v0_slider.val
    v1 = v1_slider.val
    delay = delay_slider.val

    results = nat_heu(v0, v1, delay)
    if len(results) == 0:
        abc = 1
    else:
        q1, nat_q1, makespan, nat_t, approx_ratio, qx, qy = results
        scatter_q1.set_offsets(q1)
        scatter_p1.set_offsets(p1)
        scatter_nat_q1.set_offsets(nat_q1)
        plot_handoff.set_xdata(qx)
        plot_handoff.set_ydata(qy)
        plot_opt.set_xdata([p1[0], q1[0], t[0]])
        plot_opt.set_ydata([p1[1], q1[1], t[1]])
        plot_nat.set_xdata([p1[0], nat_q1[0], t[0]])
        plot_nat.set_ydata([p1[1], nat_q1[1], t[1]])
        ax.set_title(r'Approx-ratio: ' + ('%.3f' % approx_ratio) + r', makespan: ' + ('%.3f' % makespan) +
                 r', natural: ' + ('%.3f' % nat_t))
    #l.set_ydata(amp*np.sin(2*np.pi*freq*t))

    fig.canvas.draw_idle()


def reset(event):
    ax_v0.reset()
    ax_v1.reset()


def on_press(event):
    if event.inaxes != ax: return
    #artist = event.artist
    x, y = event.xdata, event.ydata
    p1[0] = x
    p1[1] = y
    v0 = v0_slider.val
    v1 = v1_slider.val
    delay = delay_slider.val

    results = nat_heu(v0, v1, delay)
    if len(results) == 0:
        abc = 1
    else:
        q1, nat_q1, makespan, nat_t, approx_ratio, qx, qy = results
        scatter_p1.set_offsets(p1)
        scatter_q1.set_offsets(q1)
        scatter_nat_q1.set_offsets(nat_q1)
        plot_handoff.set_xdata(qx)
        plot_handoff.set_ydata(qy)
        plot_opt.set_xdata([p1[0], q1[0], t[0]])
        plot_opt.set_ydata([p1[1], q1[1], t[1]])
        plot_nat.set_xdata([p1[0], nat_q1[0], t[0]])
        plot_nat.set_ydata([p1[1], nat_q1[1], t[1]])
        ax.set_title(r'Approx-ratio: ' + ('%.3f' % approx_ratio) + r', makespan: ' + ('%.3f' % makespan) +
                     r', natural: ' + ('%.3f' % nat_t))
    fig.canvas.draw_idle()


def colorfunc(label):
    scatter_s.set_color(label)
    #plot1.set_color(label)
    #plot1.set_color(label)
    fig.canvas.draw_idle()


if __name__ == '__main__':
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.axis('equal')
    #ax.axis([-2,2,-10,10])
    # Create data for the plot
    s = np.asarray([0., 0.])
    t = np.asarray([1., 0.])
    p1 = np.asarray([1., 1.])
    #delta_f = 5.0
    #s = a0 * np.sin(2 * np.pi * f0 * t)
    v0_init = 1
    v1_init = 2
    delay_init = 0.4

    q1, nat_q1, makespan, nat_t, approx_ratio, qx, qy = nat_heu(v0_init, v1_init, delay_init)

    # Plot handler
    #plt.scatter([-1,1,-1,1], [-1,-1,1.5,1.5], c=[[1, 1, 1]])
    scatter_s = plt.scatter(s[0], s[1], c=[[1,0,0]])
    scatter_t = plt.scatter(t[0], t[1], c=[[0, 1, 0]])
    scatter_p1 = plt.scatter(p1[0], p1[1], c=[[0, 0, 0]])
    scatter_q1 = plt.scatter(q1[0], q1[1], c=[[1, 0.7, 0.2]])
    scatter_nat_q1 = plt.scatter(nat_q1[0], nat_q1[1], c=[[1, 0.2, 0.7]])
    plot_handoff, = plt.plot(qx, qy, c=[1, 0.7, 0.2])
    plot_opt, = plt.plot([p1[0], q1[0], t[0]], [p1[1], q1[1], t[1]], c=[1, 0.9, 0.5], zorder=0, linestyle='--')
    plot_nat, = plt.plot([p1[0], nat_q1[0], t[0]], [p1[1], nat_q1[1], t[1]], c=[1, 0.5, 0.9], zorder=0, linestyle='--')
    #ax.xaxis.set_label_coords(0.5, -1)
    ax.set_title(r'Approx-ratio: ' + ('%.3f' % approx_ratio) + r', makespan: ' + ('%.3f' % makespan) +
                  r', natural: ' + ('%.3f' % nat_t))
    ax.margins(x=0)

    axcolor = 'lightgoldenrodyellow'
    #axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    #axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    ax_v0 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
    ax_v1 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    ax_delay = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

    v0_slider = Slider(ax_v0, 'v0', 0.1, 30.0, valinit=v0_init)
    v1_slider = Slider(ax_v1, 'v1', 0.1, 10.0, valinit=v1_init)
    delay_slider = Slider(ax_delay, 'delay', 0, 10, valinit=delay_init)

    v0_slider.on_changed(update)
    v1_slider.on_changed(update)
    delay_slider.on_changed(update)

    resetax = plt.axes([0.8, 0.00, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    button.on_clicked(reset)

    rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
    radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)

    radio.on_clicked(colorfunc)
    fig.canvas.mpl_connect('button_press_event', on_press)
    plt.show()