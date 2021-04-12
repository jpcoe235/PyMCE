def calc_ekin_tr():
    dr_com = np.zeros(3)
    v_com = np.zeros_like(dr_com)

    for i in range(geo.natoms):
        for j in range(3):
            dr_com[j] += np.sqrt(geo.masses[i]) * (q_new_ref[j, i] - q_old_ref[j, i]) / totmass

    for i in range(3):
        v_com[j] = dr_com[j] / dt

    # v_com = v_com / totmass

    idf = 0
    for i in range(geo.natoms):
        for j in range(3):
            q[idf] = q_new_ref[j, i] - dr_com[j] * np.sqrt(geo.masses[i])
            p[idf] = p_new_ref[j, i] - v_com[j] * np.sqrt(geo.masses[i])
            idf += 1

    for i in range(3):
        ekin_tr = ekin_tr + 0.5 * (np.sqrt(totmass) * v_com[i]) ** 2.0