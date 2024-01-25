import numpy as np
from initialize_traj import trajectory


class full_trajectory():
    def __init__(self, end_time, dt, ndim, nstates):
        npart = int(ndim / 3)
        self.n_time_p = int(end_time / dt) + 10
        print('number_of_time_points', self.n_time_p)
        self.time = np.zeros(self.n_time_p)
        self.mom = np.zeros((ndim, self.n_time_p))
        self.pos = np.zeros((ndim, self.n_time_p))
        self.amps = np.zeros((nstates, self.n_time_p), dtype=np.complex128)
        self.derivs = np.zeros((ndim, nstates, nstates, self.n_time_p))
        self.HE = np.zeros((nstates, nstates, self.n_time_p), dtype=np.complex128)
        self.widths = np.zeros((ndim, self.n_time_p))
        self.phases = np.zeros(self.n_time_p)
        self.nstates = nstates
        self.ndim = ndim
        self.epot = np.zeros((nstates, self.n_time_p))
        self.mass = np.zeros(ndim)

    def get_full_phase(self):
        return self.phases

    def get_ndim(self):
        return self.ndim

    def get_widths(self):
        return self.widths

    def get_widths_time(self, time):
        return self.widths[:, time]

    def get_time_phases(self, NN):
        return self.phases[NN]

    def get_full_mom(self):
        return self.mom

    def get_full_pos(self):
        return self.pos

    def get_full_amps(self):
        return self.amps

    def get_full_derivs(self):
        return self.derivs

    def get_full_time(self):
        return self.time

    def get_mom_time(self, N):
        return self.mom[:, N]

    def get_pos_time(self, N):
        return self.pos[:, N]

    def get_derivs_time(self, N):
        return self.derivs[:, :, :, N]

    def get_derivs_time_force(self, i, N):
        return self.derivs[:, i, i, N]

    def get_derivs_time_states(self, i, j, N):
        return self.derivs[:, i, j, N]

    def get_mass(self):
        return self.mass

    def get_amps_time(self, N):
        return self.amps[:, N]

    def get_amps(self):
        return self.amps

    def get_amps_time_state(self, i, N):
        return self.amps[i, N]

    def get_HE(self):
        return self.HE

    def get_HE_time(self, N):
        return self.HE[:, :, N]

    def get_potential_full(self):
        return self.epot

    def get_potential_traj(self, i, N):
        return self.epot[i, N]

    def set_mass(self, value):
        self.mass = value

    def set_HE(self, HE):
        self.HE = HE

    def set_potential_full(self, value):
        self.epot = value

    def set_potential_time(self, value, N):
        self.epot[:, N] = value

    def set_full_phases(self, phase):
        self.phases = phase

    def set_widths(self, width):
        for i in range(self.n_time_p):
            self.widths[:, i] = width

    def set_time_phase(self, phase, N):
        self.phases[N] = phase

    def set_full_mom(self, mom):
        self.mom = mom

    def set_full_pos(self, pos):
        self.pos = pos

    def set_full_amps(self, amps):
        self.amps = amps

    def set_full_derivs(self, derivs):
        self.derivs = derivs

    def set_full_time(self, time):
        self.time = time

    def set_time_time(self, N, time1):
        self.time[N] = time1

    def set_mom_time(self, N, mom):
        self.mom[:, N] = mom

    def set_pos_time(self, N, pos):
        self.pos[:, N] = pos

    def set_derivs_time(self, N, derivs):
        self.derivs[:, :, :, N] = derivs

    def set_amps_time(self, N, amps):
        self.amps[:, N] = amps

    def set_HE_time(self, N, HE):
        self.HE[:, :, N] = HE

    def set_update_time(self, N, time, T):
        self.pos[:, N] = T.getposition_traj()
        self.mom[:, N] = T.getmomentum_traj()
        self.amps[:, N] = T.getamplitude_traj()
        derivs = np.zeros((T.ndim, T.nstates, T.nstates))
        for i in range(T.nstates):
            derivs[:, i, i] = T.getforce_traj(i)
            for j in range(T.nstates):
                derivs[:, i, j] = T.getcoupling_traj(i, j)

        self.derivs[:, :, :, N] = derivs
        self.HE[:, :, N] = T.getHE_traj()
        self.time[N] = time
        self.phases[N] = T.getphase_traj()

    def get_traj_force(self, N):
        nst = self.nstates
        f1 = np.zeros(self.ndim)
        f2 = np.zeros(self.ndim)

        E = np.zeros(nst)
        a = np.zeros(nst, dtype=np.complex128)
        for i in range(nst):
            E[i] = self.get_potential_traj(i, N)
            a[i] = self.get_amps_time_state(i, N)

        for i in range(nst):
            f1 += self.get_derivs_time_force(i, N) * np.abs(a[i]) ** 2
        for i in range(nst):
            for j in range(i + 1, nst):
                tmp = 2.0 * np.real(np.conj(a[i]) * a[j]) * (E[i] - E[j])
                f2 += tmp * self.get_derivs_time_states(i, j, N)

        fvec = f1 + f2

        return fvec

    def get_calc_HE_traj(self, N):
        nstates = self.nstates
        HE = np.zeros((nstates, nstates), dtype=np.complex128)
        for i in range(nstates):
            HE[i, i] = self.get_potential_traj(i, N)
            for j in range(i + 1, nstates):
                HE[i, j] = -1j * self.coupdotvel(i, j, N)
                HE[j, i] = -HE[i, j]
        self.set_HE_time(N, HE)
        return self.HE[:, :, N]

    def getvelocity_time(self, N):
        V = self.get_mom_time(N) / self.get_mass()
        return V

    def coupdotvel(self, i, j, N):
        if i == j:
            coup = 0.0
        else:
            coup = np.dot(self.getvelocity_time(N), self.get_derivs_time_states(i, j, N))

        return coup
