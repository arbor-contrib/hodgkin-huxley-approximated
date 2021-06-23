#!/usr/bin/env python3
"""Hodgkin-Huxley neuron simulation with approximations for gating variable steady-states and time constants

Follows exercise 4, chapter 2 of Eugene M. Izhikevich: Dynamical Systems in Neuroscience

Sebastian Schmitt, 2021
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import arbor


class HodgkinHuxleyApproximated(arbor.recipe):
    def __init__(self, catalogue, probes):
        """Initialize the recipe
        catalogue -- catalogue of mechanisms
        probes -- list of probes
        """

        # (4.1) The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)

        self.the_probes = probes
        self.the_props = arbor.neuron_cable_properties()
        self.the_cat = catalogue
        self.the_props.register(self.the_cat)

    def num_cells(self):
        return 1

    def cell_kind(self, gid):
        return arbor.cell_kind.cable

    def cell_description(self, gid):

        tree = arbor.segment_tree()

        radius = 0.1
        tree.append(arbor.mnpos, arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius), tag=1)

        labels = arbor.label_dict({'soma':   '(tag 1)',
                                   'midpoint': '(location 0 0.5)'})

        decor = arbor.decor()
        decor.set_property(Vm=0)

        decor.set_ion("na", rev_pot=120)
        decor.set_ion("k", rev_pot=-12)

        decor.paint('"soma"', arbor.mechanism("hh_approx",
                                              {"calc_initial": 0,
                                               "m_initial": 0.05,
                                               "n_initial": 0.32,
                                               "h_initial": 0.60})
        )

        I_stimulus = [[0, 0]]

        def add_stimulus(I, start, duration):
            decor.place('"midpoint"', arbor.iclamp(
                start, duration, I), "iclamp")

            I_stimulus.append([start, 0])
            I_stimulus.append([start, I])
            I_stimulus.append([start+duration, I])
            I_stimulus.append([start+duration, 0])

        area = 4 * np.pi * (radius * 1e-6)**2

        # convert 4 uA/cm^2 to total current in nA
        I = (4*1e-6/0.01**2 * area)/1e-9
        add_stimulus(I, 2, 0.5)

        # convert 15 uA/cm^2 to total current in nA
        I = (15e-6/0.01**2 * area)/1e-9
        add_stimulus(I, 10, 0.5)

        I_stimulus.append([20, 0])
        self.I_stimulus = np.array(I_stimulus)

        cell = arbor.cable_cell(tree, labels, decor)

        return cell

    def probes(self, gid):
        return self.the_probes

    def global_properties(self, kind):
        return self.the_props


def plot_membrane_voltage(ax, t, v):
    """Plot simulation result: membrane potential.

    ax -- matplotlib axes to be plotted on
    t -- list of simulation times
    v -- list of membrane values at simulation times
    """

    ax.plot(t, v, label='membrane voltage')
    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel('$v$ (mV)')
    ax.axhline(0, linestyle='dashed')
    ax.legend()


def plot_gating_variable_activations(ax, t, m, n, h):
    """Plot simulation result: gating variables.

    ax -- matplotlib axes to be plotted on
    t -- list of simulation times
    m -- list of gating variable m activations
    n -- list of gating variable n activations
    h -- list of gating variable h activations
    """

    ax.plot(t, m, label='$m$')
    ax.plot(t, n, label='$n$')
    ax.plot(t, h, label='$h$')
    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel('activation')
    ax.legend()


def plot_conductances(ax, t, g_na, g_k):
    """Plot simulation result: conductances.

    ax -- matplotlib axes to be plotted on
    t -- list of simulation times
    g_na -- list of Na channel conductances
    g_k -- list of K channel conductances
    """

    ax.plot(t, g_k, label=r'$g_\mathregular{K}$')
    ax.plot(t, g_na, label=r'$g_\mathregular{Na}$')

    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel('$g$ (S/cm$^2$)')
    ax.legend()


def plot_currents(ax, t, I_na, I_k, I_total):
    """Plot simulation result: currents.

    ax -- matplotlib axes to be plotted on
    t -- list of simulation times
    I_na -- list of Na channel currents
    I_k -- list of K channel currents
    I_total -- sum of membrane currents
    """

    ax.plot(t, I_k, label=r'$I_\mathregular{K}$')

    ax.plot(t, I_na, label=r'$I_\mathregular{Na}$')

    ax.plot(t, I_total,
            label=r'$I_\mathregular{Na} + I_\mathregular{K} + I_\mathregular{L}$')

    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel(r'I ($\mu$A/cm$^2$)')
    ax.legend()


def plot_current_stimulus(ax, I):
    """Plot simulation result: external current stimulus.

    ax -- matplotlib axes to be plotted on
    I -- numpy array of pairs of time and current values of external stimulus
    """

    ax.plot(I[:, 0], I[:, 1], label=r'$I_\mathregular{ext}$')

    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel('I (nA)')
    ax.legend()


def plot_gating_variable_time_constants(ax, t, tau_m, tau_n, tau_h):
    """Plot simulation result: gating variable time constants.

    ax -- matplotlib axes to be plotted on
    t -- list of simulation times
    m -- list of gating variable m time constants
    n -- list of gating variable n time constants
    h -- list of gating variable h time constants
    """

    ax.plot(t, tau_m, label=r'$\tau_m$')
    ax.plot(t, tau_n, label=r'$\tau_n$')
    ax.plot(t, tau_h, label=r'$\tau_h$')

    ax.set_xlabel('$t$ (ms)')
    ax.set_ylabel(r'$\tau$ (ms)')
    ax.legend()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Hodgkin-Huxley with approximations for gating variable steady-states and time constants')

    parser.add_argument(
        '--catalogue', help="name of catalogue file library", default="hh_approx-catalogue.so")
    parser.add_argument('--show', help="show plot",
                        action="store_true", default=False)
    parser.add_argument('--save', help="save to given file name")

    # parse the command line arguments
    args = parser.parse_args()

    if not args.show and not args.save:
        print("Neither --show nor --save selected, "
              "simulation will run but no output will be produced.")

    # load custom catalogue
    catalogue = arbor.load_catalogue(args.catalogue)

    # set up probes
    probes = [arbor.cable_probe_membrane_voltage('"midpoint"'),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "m"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "n"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "h"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "tau_m"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "tau_n"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "tau_h"),
              arbor.cable_probe_ion_current_density('"midpoint"', "na"),
              arbor.cable_probe_ion_current_density('"midpoint"', "k"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "gna"),
              arbor.cable_probe_density_state('"midpoint"', "hh_approx", "gk"),
              arbor.cable_probe_total_ion_current_density('"midpoint"')
              ]

    # instantiate recipe
    recipe = HodgkinHuxleyApproximated(catalogue, probes)

    # create a default execution context and a default domain decomposition
    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)

    # configure the simulation and handles for the probes
    sim = arbor.simulation(recipe, domains, context)

    # time step for simulation and sampling in ms
    dt = 0.01

    membrane_handle = sim.sample((0, 0), arbor.regular_schedule(dt))
    m_handle = sim.sample((0, 1), arbor.regular_schedule(dt))
    n_handle = sim.sample((0, 2), arbor.regular_schedule(dt))
    h_handle = sim.sample((0, 3), arbor.regular_schedule(dt))

    tau_m_handle = sim.sample((0, 4), arbor.regular_schedule(dt))
    tau_n_handle = sim.sample((0, 5), arbor.regular_schedule(dt))
    tau_h_handle = sim.sample((0, 6), arbor.regular_schedule(dt))

    I_na_handle = sim.sample((0, 7), arbor.regular_schedule(dt))
    I_k_handle = sim.sample((0, 8), arbor.regular_schedule(dt))

    g_na_handle = sim.sample((0, 9), arbor.regular_schedule(dt))
    g_k_handle = sim.sample((0, 10), arbor.regular_schedule(dt))

    I_total_handle = sim.sample((0, 11), arbor.regular_schedule(dt))

    # run the simulation for 20 ms
    sim.run(tfinal=20, dt=dt)

    # sample results
    t = sim.samples(membrane_handle)[0][0][:, 0]
    v = sim.samples(membrane_handle)[0][0][:, 1]

    m = sim.samples(m_handle)[0][0][:, 1]
    n = sim.samples(n_handle)[0][0][:, 1]
    h = sim.samples(h_handle)[0][0][:, 1]

    tau_m = sim.samples(tau_m_handle)[0][0][:, 1]
    tau_n = sim.samples(tau_n_handle)[0][0][:, 1]
    tau_h = sim.samples(tau_h_handle)[0][0][:, 1]

    I_na = sim.samples(I_na_handle)[0][0][:, 1]
    I_k = sim.samples(I_k_handle)[0][0][:, 1]

    g_na = sim.samples(g_na_handle)[0][0][:, 1]
    g_k = sim.samples(g_k_handle)[0][0][:, 1]

    I_total = sim.samples(I_total_handle)[0][0][:, 1]

    # plot
    linestyle_cycler = cycler('linestyle', ['-', '--', ':', '-.'])
    plt.rc('axes', prop_cycle=linestyle_cycler)

    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    ax0, ax1, ax2, ax3, ax4, ax5 = fig.subplots(6)

    plot_membrane_voltage(ax0, t, v)
    plot_gating_variable_activations(ax1, t, m, n, h)
    plot_conductances(ax2, t, g_na, g_k)
    plot_currents(ax3, t, I_na, I_k, I_total)
    plot_current_stimulus(ax4, recipe.I_stimulus)
    plot_gating_variable_time_constants(ax5, t, tau_m, tau_n, tau_h)

    if args.save:
        fig.savefig(args.save)

    if args.show:
        plt.show()
