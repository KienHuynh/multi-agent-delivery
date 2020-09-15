import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial


def point_comment(ax, label, x, y):
    ax.annotate(label,
                xy=(x[-1], y[-1]),
                xytext=(1.02 * x[-1], y[-1]),
                )


def solo_approx(xy_, v0_, v1_):
    """
    Find the setting with the worst approximation result when no hand-off is allowed
    :param xy_
    :param v0_:
    :param v1_:
    :return:
    """
    figure_path = './fig/solo/'

    p0 = np.asarray([0, 0])
    p1 = np.asarray([1, 1])
    v0 = 1
    v1 = 2
    s = np.asarray([0, 0])
    t = np.asarray([1, 0])

    # for pairwise computation
    eps = 0.1
    # xyp = np.repeat(xy_,)

    d0 = np.sqrt(np.sum((xy_ - p0) ** 2, axis=1))
    d0 = np.reshape(d0, (xy_.shape[0], 1))
    d0 = np.repeat(d0, xy_.shape[0], axis=1)
    d1 = scipy.spatial.distance.cdist(xy_, xy_)
    a1_s_t_arr = np.sqrt(np.sum((xy_ - s) ** 2, axis=1)) + t[0]
    # triu = np.triu_indices(xy_.shape[0])
    # d1[triu] = -1

    # solo_dist = scipy.spatial.distance.cdist(xy_, s) + t
    plt.ion()
    # fig = plt.figure(1)

    # A structure with the following data
    # 0: v0
    # 1: v1
    # 2: position of p1 in xy_
    # 3: position of the best handoff point in xy_
    # 4: approx ratio = (6) : (5)
    # 5: best makespan
    # 6: solo makespan
    hepta_list = np.zeros((int(v0_.shape[0] * (v1_.shape[0] - 1) / 2), 7))

    hepta_count = 0
    fig, ax = plt.subplots()
    for v0 in v0_:
        for v1 in v1_:
            print([v0, v1])
            if v1 <= v0:
                continue

            # for p1 in xy_:

            # Compute n^2 pair-wise distance
            eq = np.argwhere(np.bitwise_and(
                np.abs(d1 / v1 - d0 / v0) < 0.01, d0 < 1
            ))

            un = np.unique(eq[:, 1])
            hepta = np.zeros((7,))
            hepta[0] = v0
            hepta[1] = v1

            for u in un:
                p1 = xy_[u, :]
                indices = eq[eq[:, 1] == u, 0]
                handoff_points = xy_[indices, :]

                # Check if (a1 -> handoff -> t) is actually faster than (a0 -> t)
                makespan = (
                                   np.sqrt(np.sum((handoff_points - p1) ** 2, axis=1)) +
                                   np.sqrt(np.sum((handoff_points - t) ** 2, axis=1))
                           ) / v1
                min_makespan_ind = np.argmin(makespan)
                min_makespan = makespan[min_makespan_ind]
                min_handoff_ind = indices[min_makespan_ind]
                min_handoff_point = handoff_points[min_makespan_ind]

                a0_t = np.sqrt(np.sum(t - p0) ** 2) / v0
                # handoff_points = handoff_points[makespan < a0_t, :]
                if min_makespan >= a0_t:
                    continue

                a1_s_t = a1_s_t_arr[u] / v1
                solo_t = min(a1_s_t, a0_t)

                approx_ratio = solo_t / min_makespan

                if (hepta[4] < approx_ratio):
                    hepta[2] = u
                    hepta[3] = min_handoff_ind
                    hepta[4] = approx_ratio
                    hepta[5] = min_makespan
                    hepta[6] = solo_t

            hepta_list[hepta_count, :] = hepta
            hepta_count += 1

    worst_approx = hepta_list[np.argsort(hepta_list[:, 4])[-1], :]
    print(np.argsort(hepta_list[:, 4])[-1])
    print(worst_approx)

    for i in range(hepta_list.shape[0]):
        v0 = hepta_list[i][0]
        v1 = hepta_list[i][1]
        u = int(hepta_list[i][2])
        p1 = xy_[u]
        min_handoff_point = xy_[int(hepta_list[i][3]), :]
        approx_ratio = hepta_list[i][4]
        min_makespan = hepta_list[i][5]
        solo_t = hepta_list[i][6]

        ax.grid(color=[0.8, 0.8, 0.8], zorder=0)

        s_scatter = ax.scatter(s[0], s[1], c=[[1, 0, 0]], zorder=3)
        t_scatter = ax.scatter(t[0], t[1], c=[[0, 1, 0]], zorder=3)
        p1_scatter = ax.scatter(xy_[u, 0], xy_[u, 1], c=[[0, 0, 0]], zorder=3)

        # p1 -> q1 -> t
        ax.arrow(xy_[u, 0],
                 xy_[u, 1],
                 min_handoff_point[0] - xy_[u, 0],
                 min_handoff_point[1] - xy_[u, 1],
                 width=0.005,
                 color=[0, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(min_handoff_point[0],
                 min_handoff_point[1],
                 t[0] - min_handoff_point[0],
                 t[1] - min_handoff_point[1],
                 width=0.005,
                 color=[0, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(s[0],
                 s[1],
                 min_handoff_point[0] - s[0],
                 min_handoff_point[1] - s[1],
                 width=0.005,
                 color=[0, 0, 0.5],
                 zorder=2,
                 length_includes_head=True)

        # ax.plot([t[0], min_handoff_point[0]], [t[1], min_handoff_point[1]], c=[0, 0, 1])
        # ax.scatter(xy_[:,0], xy_[:,1], c=[0.8, 0.8, 0.8], marker = '+', zorder=1)
        handoff_scatter = ax.scatter(min_handoff_point[0], min_handoff_point[1], c=[[1, 0.7, 0.2]], zorder=3)

        ax.legend((s_scatter, t_scatter, p1_scatter, handoff_scatter),
                  ('s' + str(s) + ' ' + str(v0),
                   't' + str(t),
                   'a1' + str(p1) + ' ' + str(v1),
                   'q' + str(min_handoff_point)),
                  scatterpoints=1)

        ax.axis('equal')
        ax.set_xlabel(r'Approx-ratio: ' + ('%.3f' % approx_ratio) + r', makespan: ' + ('%.3f' % min_makespan) +
                      r', solo: ' + ('%.3f' % solo_t))
        plt.gca().set_aspect("equal")
        plt.draw()
        plt.show()
        # plt.savefig(figure_path + ('%04d.pdf' % i))
        plt.cla()


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


def nat_approx(xy_, v0_, v1_):
    figure_path = './fig/nat/'

    p0 = np.asarray([0, 0])
    p1 = np.asarray([1, 1])
    v0 = 1
    v1 = 2
    s = np.asarray([0, 0])
    t = np.asarray([1, 0])

    # for pairwise computation
    eps = 0.1
    # xyp = np.repeat(xy_,)

    d0 = np.sqrt(np.sum((xy_ - p0) ** 2, axis=1))
    d0 = np.reshape(d0, (xy_.shape[0], 1))
    d0 = np.repeat(d0, xy_.shape[0], axis=1)
    d1 = scipy.spatial.distance.cdist(xy_, xy_)
    a1_s_t_arr = np.sqrt(np.sum((xy_ - s) ** 2, axis=1)) + t[0]

    # solo_dist = scipy.spatial.distance.cdist(xy_, s) + t
    plt.ion()
    # fig = plt.figure(1)

    # A structure with the following data
    # 0: v0
    # 1: v1
    # 2: position of p1 in xy_
    # 3: position of the best handoff point in xy_
    # 4: xq1
    # 5: approx ratio = (6) : (5)
    # 6: best makespan
    # 7: natural heuristic makespan
    hepta_list = np.zeros((int(v0_.shape[0] * (v1_.shape[0] - 1) / 2), 8))

    hepta_count = 0
    fig, ax = plt.subplots()
    for v0 in v0_:
        for v1 in v1_:
            print([v0, v1])
            if v1 <= v0:
                continue

            # for p1 in xy_:

            # Compute n^2 pair-wise distance
            bin_eq = np.bitwise_and(np.abs(d1 / v1 - d0 / v0) < 0.01, d0 < 1)
            eq = np.argwhere(bin_eq)

            un = np.unique(eq[:, 1])
            hepta = np.zeros((8,))
            hepta[0] = v0
            hepta[1] = v1

            for u in un:
                p1 = xy_[u, :]
                indices = eq[eq[:, 1] == u, 0]
                handoff_points = xy_[indices, :]

                # Check if (a1 -> handoff -> t) is actually faster than (a0 -> t)
                makespan = (
                                   np.sqrt(np.sum((handoff_points - p1) ** 2, axis=1)) +
                                   np.sqrt(np.sum((handoff_points - t) ** 2, axis=1))
                           ) / v1
                min_makespan_ind = np.argmin(makespan)
                min_makespan = makespan[min_makespan_ind]
                min_handoff_ind = indices[min_makespan_ind]
                min_handoff_point = handoff_points[min_makespan_ind]

                # natural heuristic computation
                xq1 = solve_quad(p1[0], p1[1], s[0], t[0], v0, v1)
                if (xq1 == -1):
                    continue
                nat_q1 = np.asarray([xq1, 0])
                nat_t = (
                                np.sqrt(np.sum((nat_q1 - p1) ** 2)) +
                                np.sqrt(np.sum((nat_q1 - t) ** 2))
                        ) / v1

                a0_t = np.sqrt(np.sum(t - p0) ** 2) / v0

                if min_makespan >= a0_t:
                    continue

                a1_s_t = a1_s_t_arr[u] / v1

                approx_ratio = nat_t / min_makespan

                if (hepta[4] < approx_ratio):
                    hepta[2] = u
                    hepta[3] = min_handoff_ind
                    hepta[4] = xq1
                    hepta[5] = approx_ratio
                    hepta[6] = min_makespan
                    hepta[7] = nat_t

            hepta_list[hepta_count, :] = hepta
            hepta_count += 1

    worst_approx = hepta_list[np.argsort(hepta_list[:, 5])[-1], :]
    print(np.argsort(hepta_list[:, 5])[-1])
    print(worst_approx)

    for i in range(hepta_list.shape[0]):
        v0 = hepta_list[i][0]
        v1 = hepta_list[i][1]
        u = int(hepta_list[i][2])
        p1 = xy_[u]
        min_handoff_point = xy_[int(hepta_list[i][3]), :]
        q1 = np.asarray([hepta_list[i][4], 0])
        approx_ratio = hepta_list[i][5]
        min_makespan = hepta_list[i][6]
        nat_t = hepta_list[i][7]

        ax.grid(color=[0.8, 0.8, 0.8], zorder=0)

        s_scatter = ax.scatter(s[0], s[1], c=[[1, 0, 0]], zorder=3)
        t_scatter = ax.scatter(t[0], t[1], c=[[0, 1, 0]], zorder=3)
        p1_scatter = ax.scatter(xy_[u, 0], xy_[u, 1], c=[[0, 0, 0]], zorder=3)

        # OPT
        # p1 -> q1 -> t
        ax.arrow(xy_[u, 0],
                 xy_[u, 1],
                 min_handoff_point[0] - xy_[u, 0],
                 min_handoff_point[1] - xy_[u, 1],
                 width=0.005,
                 color=[0, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(min_handoff_point[0],
                 min_handoff_point[1],
                 t[0] - min_handoff_point[0],
                 t[1] - min_handoff_point[1],
                 width=0.005,
                 color=[0, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(s[0],
                 s[1],
                 min_handoff_point[0] - s[0],
                 min_handoff_point[1] - s[1],
                 width=0.005,
                 color=[0, 0, 0.5],
                 zorder=2,
                 length_includes_head=True)

        # Natural heuristic
        ax.arrow(xy_[u, 0],
                 xy_[u, 1],
                 q1[0] - xy_[u, 0],
                 q1[1] - xy_[u, 1],
                 width=0.005,
                 color=[1, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(q1[0],
                 q1[1],
                 t[0] - q1[0],
                 t[1] - q1[1],
                 width=0.005,
                 color=[1, 0, 1],
                 zorder=2,
                 length_includes_head=True)

        ax.arrow(s[0],
                 s[1],
                 q1[0] - s[0],
                 q1[1] - s[1],
                 width=0.005,
                 color=[0.5, 0, 0.5],
                 zorder=2,
                 length_includes_head=True)

        handoff_scatter = ax.scatter(min_handoff_point[0], min_handoff_point[1], c=[[1, 0.7, 0.2]], zorder=3)
        q1_scatter = ax.scatter(q1[0], q1[1], c=[[1, 0.2, 0.7]], zorder=3)

        ax.legend((s_scatter, t_scatter, p1_scatter, handoff_scatter, q1_scatter),
                  ('s' + str(s) + ' ' + str(v0),
                   't' + str(t),
                   'a1' + str(p1) + ' ' + str(v1),
                   'q' + str(min_handoff_point),
                   'nat q' + str(q1)),
                  scatterpoints=1)

        ax.axis('equal')
        ax.set_xlabel(r'Approx-ratio: ' + ('%.3f' % approx_ratio) + r', makespan: ' + ('%.3f' % min_makespan) +
                      r', natural: ' + ('%.3f' % nat_t))
        plt.gca().set_aspect("equal")
        plt.draw()
        plt.show()
        plt.savefig(figure_path + ('%04d.pdf' % i))
        plt.cla()


if __name__ == '__main__':
    np.set_printoptions(precision=3)

    x_ = np.linspace(-21, 1, 150)
    y_ = np.linspace(0, 5, 100)
    xy_ = np.meshgrid(x_, y_)
    xy_ = np.reshape(xy_, (2, 150 * 100)).T
    v0_ = np.arange(0.1, 10, 1)
    v1_ = np.arange(0.1, 10, 1)
    # solo_approx(xy_, v0_, v1_)
    nat_approx(xy_, v0_, v1_)
