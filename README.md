# Hodgkin-Huxley Approximated

[![Hodgkin-Huxley Approximated](https://github.com/schmitts/hodgkin-huxley-approximated/actions/workflows/basic.yml/badge.svg)](https://github.com/schmitts/hodgkin-huxley-approximated/actions/workflows/basic.yml)

Hodgkin-Huxley neuron simulation with approximations for gating variable steady-states and time constants.

Follows exercise 4, chapter 2 of Eugene M. Izhikevich: Dynamical Systems in Neuroscience

# Howto

Install Arbor from source and make sure `build-catalogue` from the Arbor scripts directory
is in `PATH` or modify the `Makefile` accordingly.

```shell
make
```

![Hodgkin-Huxley Approximated](hh_approximated.svg)