import numpy as np
import json
from itertools import permutations
import matplotlib.pyplot as plt
import copy
import time


def gen(min_x, max_x, min_y, max_y, min_speed, max_speed, m, k):
    scenario = {}
    agents = []
    packages = []

    scenario['k'] = k
    scenario['m'] = m
    for mm in range(m):
        package_mm = {}
        x = np.random.uniform(min_x, max_x)
        y = np.random.uniform(min_y, max_y)
        xy = [x, y]
        package_mm['s'] = xy
        x = np.random.uniform(min_x, max_x)
        y = np.random.uniform(min_y, max_y)
        xy = [x, y]
        package_mm['t'] = xy
        packages.append(package_mm)

    for kk in range(k):
        agent_kk = {}
        x = np.random.uniform(min_x, max_x)
        y = np.random.uniform(min_y, max_y)
        v = np.random.uniform(min_speed, max_speed)
        xy = [x, y]
        agent_kk['loc'] = xy
        agent_kk['v'] = v
        agents.append(agent_kk)

    agents = sorted(agents, key=lambda i: i['v'])

    scenario['agents'] = agents
    scenario['packages'] = packages

    return scenario


def write_input_file(scenario, save_path):
    txt_save_path = save_path + '.txt'
    json_save_path = save_path + '.json'

    with open(txt_save_path, 'w') as f:
        f.write('TWODIM EUCLID DISCRETE' + '\n')

        f.write('SINGLE_POINT' + '\n')
        f.write(str(scenario['m']) + '\n')
        for mm in range(scenario['m']):
            package = scenario['packages'][mm]
            f.write(str(package['s'][0]) + ' ' + str(package['s'][1]) + ' ' + str(mm) + '\n')

        f.write('SINGLE_POINT' + '\n')
        f.write(str(scenario['m']) + '\n')
        for mm in range(scenario['m']):
            package = scenario['packages'][mm]
            f.write(str(package['t'][0]) + ' ' + str(package['t'][1]) + ' ' + str(mm) + '\n')

        f.write(str(scenario['k']) + '\n')
        for kk in range(scenario['k']):
            agent = scenario['agents'][kk]
            f.write(str(agent['loc'][0]) + ' ' + str(agent['loc'][1]) + ' ' + str(agent['v']) + '\n')

        f.write('APGRID' + '\n')
        f.write('-20 20 16' + '\n')
        f.write('-20 20 16' + '\n')

    scenario_dict = {'scenario': scenario}
    with open(json_save_path, 'w') as f:
        f.write(json.dumps(scenario_dict))


def timing(a, b, v):
    return np.sqrt(np.sum((a - b) ** 2)) / v


def solve_quad(a, b, c):
    delt = b ** 2 - 4 * a * c
    x1 = (-b + np.sqrt(delt)) / (2 * a)
    x2 = (-b - np.sqrt(delt)) / (2 * a)
    return (x1, x2)


def find_intersection(i0, i1, t0_, t1_, agent):
    # Translate i0 to (0,0) and other points accordingly
    i1 -= i0

    # Rotate every points so that i1 is on the x-axis
    r = np.sqrt(np.sum(i1 ** 2))
    sin_alpha = i1[1] / r
    cos_alpha = i1[0] / r
    r_mat = np.asarray([[cos_alpha, -sin_alpha], [sin_alpha, cos_alpha]])

    i1 = i1.dot(r_mat)
    xy = (agent[1:3] - i0).dot(r_mat)
    xa = xy[0]
    ya = xy[1]
    slope = (t1_ - t0_) / r
    quad_a = (1 - (slope ** 2) * (agent[3] ** 2))
    quad_b = -2 * xa - 2 * slope * t0_ * (agent[3] ** 2)
    quad_c = (xa ** 2) + (ya ** 2) - (agent[3] ** 2) * (t0_ ** 2)
    [x1, x2] = solve_quad(quad_a, quad_b, quad_c)

    p1 = []
    if (x1 >= 0 and x1 <= r):
        if (x2 >= 0 and x2 <= r):
            if (x1 < x2):
                p1 = [x1, 0]
            else:
                p1 = [x2, 0]
        else:
            p1 = [x1, 0]
    else:
        p1 = [x2, 0]

    p1 = np.asarray(p1)
    t_test_2 = timing(np.asarray([0, 0]), p1, 4.182276295893604) + t0_
    r_mat = np.linalg.inv(r_mat)
    p1 = p1.dot(r_mat) + i0
    t = timing(p1, agent[1:3], agent[3])
    t_test = timing(i0, p1, 4.182276295893604) + t0_

    return [p1, t]


# Not binary search
def get_intersection_index(intervals, ren_points, agent):
    # Edge case
    t0 = timing(ren_points[0], agent[1:3], agent[3])
    if (t0 < intervals[0]):
        return [ren_points[0], t0, 0]

    for i in range(len(intervals) - 1):
        i0 = ren_points[i]
        i1 = ren_points[i + 1]
        t0 = timing(i0, agent[1:3], agent[3])
        t1 = timing(i1, agent[1:3], agent[3])
        t0_ = intervals[i]
        t1_ = intervals[i + 1]
        if ((t0 - t0_) * (t1 - t1_) < 0):
            inter, t0 = find_intersection(i0, i1.copy(), t0_, t1_, agent)
            return [inter, t0, i + 1]

    # Cannot reach any point on the s-t line segment in time
    return [-1, -1, -1]


# Polynomial time algo, not the best
def nat_approx(s, t, agents):
    agents = agents[np.argsort(agents[:, 3]), :]
    a0 = agents[0, :]
    makespan = 0
    intervals = [timing(s[0, 1:3], a0[1:3], a0[3]), timing(s[0, 1:3], a0[1:3], a0[3]) + timing(s[0, 1:3], t[0, 1:3], a0[3])]
    ren_points = [s[0, 1:3], t[0, 1:3]]
    agent_used = [a0[0]]

    for a in range(1, agents.shape[0]):
        agent = agents[a]
        point, time, i = get_intersection_index(intervals, ren_points, agent)
        if (i >= 0):
            del intervals[i:]
            del ren_points[i:]
            del agent_used[i:]
            intervals += [time, timing(point, t[0, 1:3], agent[3]) + time]
            ren_points += [point, t[0, 1:3]]
            agent_used += [agent[0]]

    return intervals, agent_used, ren_points


# Create all m choose k permuations
def create_combs(m, k):
    combs = list(permutations(range(k), m))
    for c in range(len(combs)):
        comb = combs[c]
        combs[c] = []
        for mm in range(m):
            combs[c] += [[mm, comb[mm]]]

    return combs


# agent here is a dict
def update_agent(agent, time_taken):
    move_vector = agent['to'][0] - agent['loc']
    if (np.sum(move_vector) == 0):
        abc = 1
    move_vector = agent['v'] * move_vector / np.sqrt(np.sum(move_vector ** 2))
    move_vector = move_vector * time_taken
    agent['loc'] += move_vector


def update_src(src_table, agent_table):
    for i in range(src_table.shape[0]):
        mm = src_table[i, 0]
        agents = agent_table[agent_table[:, 4] == mm, :]
        agents = agents[np.argsort(agents[:, 3]), :]
        if src_table[i, -1] == 1 and agents.size > 0:
            src_table[i, 1:3] = agents[0, 1:3]


def nat_approx_m(src_table, tgt_table, agent_table, in_prog_pkg):
    agents = []
    for mm in in_prog_pkg:
        agent_table_mm = agent_table[agent_table[:, 4] == mm, :]
        if (len(agent_table_mm.shape) == 1):
            agent_table_mm.reshape(1, agent_table_mm.size)

        # Sort them by speed
        agent_table_mm = agent_table_mm[np.argsort(agent_table_mm[:, 3]), :]

        intervals, agents_used, ren_points = nat_approx(src_table[src_table[:, 0] == mm, :],
                                              tgt_table[src_table[:, 0] == mm, :],
                                              agent_table_mm)
        for kk in range(len(intervals) - 1):
            # mm: the src-target pair that this agent is assigned to
            agent = {'id': agents_used[kk],
                     'loc': agent_table_mm[agent_table_mm[:, 0] == agents_used[kk], 1:3][0],
                     'v': agent_table_mm[agent_table_mm[:, 0] == agents_used[kk], 3][0],
                     'mm': agent_table_mm[agent_table_mm[:, 0] == agents_used[kk], 4][0],
                     'to': [ren_points[kk], ren_points[kk + 1]],
                     'ts': [intervals[kk], intervals[kk + 1]]}

            agents.append(agent)

    return agents


def update_agent_table(agent_table, agents):
    for agent in agents:
        id = agent['id']
        agent_table[agent_table[:, 0] == id, 1:3] = agent['loc']
        agent_table[agent_table[:, 0] == id, 4] = agent['mm']


def live_approx(src_table, tgt_table, agent_table, m, ani1, ani2):
    ren_points = []
    current_time = 0

    src_table = np.concatenate((src_table, np.zeros((m, 1))), axis = 1)

    in_prog_pkg = list(range(m))
    agents = nat_approx_m(src_table, tgt_table, agent_table, in_prog_pkg)
    ani1 += copy.deepcopy(agents)
    finished = False
    counter = -1
    while True:
        counter += 1
        agents = sorted(agents, key=lambda i: i['ts'][0])

        current_time_ = current_time
        # Add the first travelling section of agents[0] to the animation event
        agents_0_ani = copy.deepcopy(agents[0])
        del agents_0_ani['ts'][0:]
        del agents_0_ani['to'][0:]
        agents_0_ani['ts'].append(current_time_)
        agents_0_ani['ts'].append(agents[0]['ts'][0])
        agents_0_ani['to'].append(np.copy(agents[0]['loc']))
        agents_0_ani['to'].append(np.copy(agents[0]['to'][0]))
        ani2.append(agents_0_ani)

        # Update the first agent and clear its firs task
        time_taken = agents[0]['ts'][0] - current_time
        current_time = agents[0]['ts'][0]
        #print('%.4f' % current_time, '%.4f' % time_taken)
        agents[0]['loc'] = agents[0]['to'][0]
        del agents[0]['to'][0]
        del agents[0]['ts'][0]

        # One of the package maybe picked up
        src_table[src_table[:, 0] == agents[0]['mm'], -1] = 1

        # If any of them finish both tasks, we have to rerun the nat approx
        update_state = False
        agent_0_mm_ = agents[0]['mm']
        ex_list = []
        if len(agents[0]['to']) == 0:
            update_state = True
            agents[0]['mm'] = -1

            # Remove the first task of the next agent that receives the package
            for a in range(1, len(agents)):
                if agents[a]['ts'][0] == current_time:
                    # Add animation
                    agents_0_ani = copy.deepcopy(agents[a])
                    del agents_0_ani['ts'][0:]
                    del agents_0_ani['to'][0:]
                    agents_0_ani['ts'].append(current_time_)
                    agents_0_ani['ts'].append(current_time)
                    agents_0_ani['to'].append(np.copy(agents[a]['loc']))
                    agents_0_ani['to'].append(np.copy(agents[a]['to'][0]))
                    ani2.append(agents_0_ani)
                    
                    agents[a]['loc'] = agents[a]['to'][0]
                    ex_list.append(agents[a]['id'])
                    del agents[a]['ts'][0]
                    del agents[a]['to'][0]

        # Update the location of other agents
        for a in range(1, len(agents)):
            # Already update above
            if agents[a]['id'] in ex_list:
                continue

            # Add animation
            agents_0_ani = copy.deepcopy(agents[a])
            del agents_0_ani['ts'][0:]
            del agents_0_ani['to'][0:]
            agents_0_ani['ts'].append(current_time_)
            agents_0_ani['ts'].append(current_time)
            agents_0_ani['to'].append(np.copy(agents[a]['loc']))
            update_agent(agents[a], time_taken)
            agents_0_ani['to'].append(np.copy(agents[a]['loc']))
            ani2.append(agents_0_ani)

        update_agent_table(agent_table, agents)
        update_src(src_table, agent_table)

        # Check if any package has been delivered
        mm_done = []
        for mm in in_prog_pkg:
            agents_mm = [agent for agent in agents if agent['mm'] == mm]
            pkg_delivered = True
            for agent_mm in agents_mm:
                if len(agent_mm['to']) != 0:
                    pkg_delivered = False
            if pkg_delivered:
                src_table = src_table[src_table[:, 0] != mm, :]
                tgt_table = tgt_table[tgt_table[:, 0] != mm, :]
                mm_done.append(mm)

        in_prog_pkg = [pkg for pkg in in_prog_pkg if pkg not in mm_done]

        # Every package is delivered, end the loop
        if len(in_prog_pkg) == 0:
            break

        time_table = np.zeros((len(in_prog_pkg), len(in_prog_pkg)))
        #time_table[:, 0] = in_prog_pkg
        agent_0 = agent_table[agent_table[:, 0] == agents[0]['id'], :].reshape(1, 5)
        no_contrib = True
        if update_state:
            for i in range(len(in_prog_pkg)):
                mm = in_prog_pkg[i]

                # The makespan without the newly freed agent
                agent_table_mm = agent_table[agent_table[:, 4] == mm, :]
                intervals, agents_used, _ = nat_approx(
                    src_table[src_table[:, 0] == mm, :],
                    tgt_table[src_table[:, 0] == mm, :],
                    np.copy(agent_table_mm))
                time_table[:, i] = intervals[-1]

                if mm == agent_0_mm_:
                    continue

                # The makespan with the newly freed agent
                agent_table_mm = np.concatenate((agent_0, agent_table_mm), 0)
                intervals, agents_used, _ = nat_approx(
                    src_table[src_table[:, 0] == mm, :],
                    tgt_table[src_table[:, 0] == mm, :],
                    np.copy(agent_table_mm))
                if agents[0]['id'] in agents_used:
                    #time_table[time_table[:, 0] == mm, 1] = intervals[-1]
                    time_table[i, i] = intervals[-1]
                    no_contrib = False

            if no_contrib:
                del agents[0]
                agent_table = agent_table[agent_table[:, 4] != -1, :]
            else:
                makespans = np.max(time_table, axis=1)
                best_ind = np.argmin(makespans)
                best_mm = in_prog_pkg[best_ind]
                if(best_mm == agent_0_mm_):
                    del agents[0]
                    agent_table = agent_table[agent_table[:, 4] != -1, :]
                    continue
                agent_table[agent_table[:, 0] == agents[0]['id'], 4] = best_mm
                agents = nat_approx_m(src_table, tgt_table, agent_table, in_prog_pkg)

                for a in range(len(agents)):
                    # remove all events where the agent picks up the package right where it's at
                    at_src = src_table[:, 1:3] == agents[a]['loc'].reshape(1,2)
                    at_src = np.sum(np.sum(at_src, axis=1) == 2) > 0

                    if agents[a]['ts'][0] == 0 and at_src:
                        del agents[a]['ts'][0]
                        del agents[a]['to'][0]

                    agents[a]['ts'] = [t + current_time for t in agents[a]['ts']]

    return current_time


def silly_approx(scenario, ani1, ani2):
    k = scenario['k']
    m = scenario['m']
    agents = scenario['agents']
    packages = scenario['packages']
    n_agent_p_pair = int(k / m)
    # column 0: index
    # column 1-2: (x,y)
    # column 3: v
    # column 4: boolean, check if it's used in the first phase matching yet
    agent_table = np.zeros((k, 5))

    # column 0: index
    # column 1-2: (x,y)
    src_table = np.zeros((m, 3))
    tgt_table = np.zeros((m, 3))

    for kk in range(k):
        agent_table[kk, 0] = kk
        agent_table[kk, 1:3] = agents[kk]['loc']
        agent_table[kk, 3] = agents[kk]['v']
    agent_table[:, 4] = -1

    for mm in range(m):
        src_table[mm, 0] = mm
        src_table[mm, 1:] = packages[mm]['s']
        tgt_table[mm, 0] = mm
        tgt_table[mm, 1:] = packages[mm]['t']

    agent_table_ = np.copy(agent_table)
    best_time = np.zeros((m, 1))
    best_time_ = np.zeros((m, 1))
    while agent_table[agent_table[:, 4] == -1, :].size > 0:
        # Check to see if any agent is unassigned
        agent_table_na = agent_table[agent_table[:, 4] == -1, :]
        time_table = np.zeros((agent_table_na.shape[0], m))
        time_without_new_agent = np.zeros((1, m))
        for kk in range(np.shape(agent_table_na)[0]):
            for mm in range(m):
                # Combine new one with already-assigned agents
                agent_ = agent_table[agent_table[:, 4] == mm, :]
                agent = agent_table_na[kk, :].reshape((1, 5))
                agent = np.concatenate((agent, agent_), 0)
                intervals, _, _ = nat_approx(
                    src_table[src_table[:, 0] == mm, :],
                    tgt_table[src_table[:, 0] == mm, :],
                    np.copy(agent))
                time_table[kk, mm] = intervals[-1]

        if agent_table_na.shape[0] < m:
            for mm in range(m):
                # Combine new one with already-assigned agents
                agent_ = agent_table[agent_table[:, 4] == mm, :]
                intervals, _, _ = nat_approx(
                    src_table[src_table[:, 0] == mm, :],
                    tgt_table[src_table[:, 0] == mm, :],
                    np.copy(agent_))
                time_table[0, mm] = intervals[-1]

        # Find the best assignment by brute-forcing over all combinations
        if (agent_table_na.shape[0] >= m):
            combs = create_combs(m, agent_table_na.shape[0])
            all_comb_time = np.zeros((len(combs), m))
            for c in range(len(combs)):
                for mm in range(m):
                    all_comb_time[c, mm] = time_table[combs[c][mm][1], combs[c][mm][0]]

            max_comb_time = np.max(all_comb_time, axis=1)
            best_comb_ind = np.argmin(max_comb_time)
            best_comb = combs[best_comb_ind]
            best_time = all_comb_time[best_comb_ind, :]
            if (np.sum(best_time == best_time_) != m):
                for mm in range(m):
                    if (best_time[mm] == best_time_[mm]):
                        continue
                    agent_table[agent_table_na[best_comb[mm][1], 0].astype(int), 4] = mm
                best_time_ = best_time

        else:
            combs = create_combs(agent_table_na.shape[0], m)
            all_comb_time = np.zeros((len(combs), m)) + time_without_new_agent
            for c in range(len(combs)):
                for kk in range(agent_table_na.shape[0]):
                    all_comb_time[c, combs[c][kk][1]] = time_table[combs[c][kk][0], combs[c][kk][1]]

            max_comb_time = np.max(all_comb_time, axis=1)
            best_comb_ind = np.argmin(max_comb_time)
            best_comb = combs[best_comb_ind]

            best_time = all_comb_time[best_comb_ind, :]

            # If there is no change between this loop and the previous loop
            # It means that all of the unused are agents do not help
            # We can stop the algorithm then
            no_change = True
            for c in range(len(best_comb)):
                if best_time[c] != best_time_[best_comb[c][1]]:
                    no_change = False

            if not no_change:
                for kk in range(agent_table_na.shape[0]):
                    if best_time[kk] < best_time_[best_comb[kk][1]]:
                        agent_table[agent_table_na[kk, 0].astype(int), best_comb[kk][1]] = best_comb[kk][1]
                        best_time_[best_comb[kk][1]] = best_time[kk]

        if (np.sum(agent_table_ == agent_table) == agent_table.size):
            break
        agent_table_ = np.copy(agent_table)
        abc = 1

    agent_table = agent_table[agent_table[:, 4] != -1, :]
    makespan_wo_share = np.max(best_time_)
    makespan_w_share = live_approx(np.copy(src_table), tgt_table, np.copy(agent_table), m, ani1, ani2)
    if np.abs(makespan_w_share - makespan_wo_share) > 0.001:
        abc = 1

    return makespan_wo_share, makespan_w_share


def create_line_anis(e, color):
    line_anis = []
    for i in range(len(e['to'])-1):
        p0 = e['to'][i]
        p1 = e['to'][i+1]
        t0 = e['ts'][i]
        t1 = e['ts'][i+1]
        line_anis.append({
            'p': [p0, p1],
            't': [t0, t1],
            'c': color,
            'drawn': False,
            'in_drawn': False})

    return line_anis


def events_to_anis(events):
    speeds = [s['v'] for s in events]
    max_speed = np.max(speeds)
    line_anis = []
    for i in range(len(events)):
        if (i == 28):
            abc = 1
        e = events[i]
        color = [0, 0, e['v'] / max_speed]
        line_ani = create_line_anis(e, color)
        line_anis += line_ani
    return line_anis


def draw_scenario(scenario, events, fig, ax):

    ax.axis('equal')
    for i in range(len(scenario['packages'])):
        p = scenario['packages'][i]
        s = p['s']
        t = p['t']
        ax.scatter(s[0], s[1], color='r')
        ax.scatter(t[0], t[1], color='g')
        ax.annotate('s' + str(i), (s[0] + 0.5, s[1]))
        ax.annotate('t' + str(i), (t[0] + 0.5, t[1]))

    speeds = [s['v'] for s in scenario['agents']]
    max_speed = np.max(speeds)
    for a in scenario['agents']:
        color = [0, 0, a['v'] / max_speed]
        ax.scatter(a['loc'][0], a['loc'][1], color=color)
        ax.annotate('%.2f' % a['v'], (a['loc'][0] + 0.5, a['loc'][1]))

    line_anis = events_to_anis(events)
    global_t0 = time.time()
    prev_t = 0

    max_x = ax.get_xlim()[1]
    max_y = ax.get_ylim()[1]
    time_text = ax.annotate('%.2f' % prev_t,
                            (max_x - 0.5, max_y - 2.5))

    v2 = np.sqrt(np.sum((line_anis[1]['p'][1] - line_anis[1]['p'][0]) ** 2)) / (line_anis[1]['t'][1] - line_anis[1]['t'][0])

    ani_done = False
    while not ani_done:
        curr_t = time.time() - global_t0
        time_text.set_text('%.2fs' % curr_t)

        ani_done = True
        for i in range(len(line_anis)):
            line_ani = line_anis[i]
            if line_ani['drawn']:
                continue

            if curr_t < line_ani['t'][0]:
                continue

            if not line_ani['in_drawn']:
                line_ani['in_drawn'] = True
                prev_t = line_ani['t'][0]

            ani_done = False
            next_t = curr_t
            if next_t > line_ani['t'][1]:
                next_t = line_ani['t'][1]
                line_ani['drawn'] = True

            if i == 28:
                abc = 1

            # Draw the segment from prev_t to curr_t
            line_duration = line_ani['t'][1] - line_ani['t'][0]
            old_p = ((prev_t - line_ani['t'][0]) / line_duration) * (line_ani['p'][1] - line_ani['p'][0]) + line_ani['p'][0]
            old_p = old_p.flatten()
            new_p = ((next_t - line_ani['t'][0]) / line_duration) * (line_ani['p'][1] - line_ani['p'][0]) + line_ani['p'][0]
            new_p = new_p.flatten()

            ax.plot([old_p[0], new_p[0]], [old_p[1], new_p[1]], color=line_ani['c'])

        plt.pause(0.001)

        fig.canvas.draw()
        plt.show(block=False)

        prev_t = curr_t

    plt.show(block=True)


def preprocess_events(events):
    for i in range(len(events)):
        events[i]['to'] = [events[i]['loc']] + events[i]['to']
        events[i]['ts'] = [0] + events[i]['ts']


if __name__ == '__main__':
    np.random.seed(13111991)
    save_path = './data/multi_package/%04d'
    n_ex = 20  # for each k

    min_x = -20.0
    max_x = 20.0
    min_y = -20.0
    max_y = 20.0
    min_speed = 1.0
    max_speed = 10.0
    m = 4
    k_list = [4, 6, 8]

    file_counter = 0
    # for k in k_list:
    #     for i in range(n_ex):
    #         scenario = gen(min_x, max_x, min_y, max_y, min_speed, max_speed, m, k)
    #         save_path_c = save_path % file_counter
    #         write_input_file(scenario, save_path_c)
    #         file_counter += 1

    plt.ion()
    fig, ax = plt.subplots(1, 1)
    # 47
    for ex in range(47, n_ex * len(k_list)):
        load_path_c = save_path % ex + '.json'
        scenario = {}
        with open(load_path_c, 'r') as f:
            scenario = json.load(f)
        scenario = scenario['scenario']
        events1 = []
        events2 = []
        silly_approx(scenario, events1, events2)

        preprocess_events(events1)
        #draw_scenario(scenario, events1, fig, ax)
        draw_scenario(scenario, events2, fig, ax)
        abc = 1
