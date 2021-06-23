TITLE Hodgkin-Huxley-Approximated
: Hodgkin-Huxley neuron with approximations for gating variable steady-states and time constants
: Follows exercise 4, chapter 2 of Eugene M. Izhikevich: Dynamical Systems in Neuroscience
: Sebastian Schmitt, 2021

NEURON {
	SUFFIX hh_approx
	USEION na READ ena WRITE ina
	USEION k  READ ek  WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE el
	RANGE gnabar, gkbar, gl
	RANGE v_n_half, v_m_half, v_h_half, k_n, k_m, k_h
	RANGE v_n_max, v_m_max, v_h_max, sigma_n, sigma_m, sigma_h, c_n_amp, c_m_amp, c_h_amp, c_n_base, c_m_base, c_h_base
	RANGE calc_initial, m_initial, n_initial, h_initial
}

PARAMETER {
	gnabar =   0.12   (S/cm2)
	gkbar  =   0.036  (S/cm2)
	gl     =   0.0003 (S/cm2)
	el     =  10.6    (mV)

    : Boltzmann function parameters
    v_n_half = 12 (mV)
    v_m_half = 25 (mV)
    v_h_half = 3 (mV)

    k_n = 15 (mV)
    k_m = 9 (mV)
    k_h = -7 (mV)

    : Gaussian function parameters
    v_n_max = -14 (mV)
    v_m_max = 27 (mV)
    v_h_max = -2 (mV)

    sigma_n = 50 (mV)
    sigma_m = 30 (mV)
    sigma_h = 20 (mV)

    c_n_amp = 4.7 (ms)
    c_m_amp = 0.46 (ms)
    c_h_amp = 7.4 (ms)

    c_n_base = 1.1 (ms)
    c_m_base = 0.04 (ms)
    c_h_base = 1.2 (ms)

    calc_initial = 1
	m_initial = 0.5
	n_initial = 0.5
	h_initial = 0.5
}

STATE { m h n m_inf h_inf n_inf tau_m tau_h tau_n gna gk}

BREAKPOINT {
	SOLVE states METHOD cnexp

	tau_m = tau(v, c_m_base, c_m_amp, v_m_max, sigma_m)
	tau_n = tau(v, c_n_base, c_n_amp, v_n_max, sigma_n)
	tau_h = tau(v, c_h_base, c_h_amp, v_h_max, sigma_h)

	m_inf = xinf(v, v_m_half, k_m)
	n_inf = xinf(v, v_n_half, k_n)
	h_inf = xinf(v, v_h_half, k_h)

	gk = gkbar * n^4
    gna = gnabar * m^3 * h
	ina = gna*(v - ena)

	ik  = gk*(v - ek)
	il  = gl*(v - el)
}

INITIAL {

	if (calc_initial != 0) {
	   m = xinf(v, v_m_half, k_m)
	   n = xinf(v, v_n_half, k_n)
	   h = xinf(v, v_h_half, k_h)
    } else {
	   m = m_initial
	   n = n_initial
	   h = h_initial
    }

	m_inf = m
	n_inf = n
	h_inf = h

	tau_m = tau(v, c_m_base, c_m_amp, v_m_max, sigma_m)
	tau_n = tau(v, c_n_base, c_n_amp, v_n_max, sigma_n)
	tau_h = tau(v, c_h_base, c_h_amp, v_h_max, sigma_h)
}

DERIVATIVE states {
	m' = (m_inf - m)/tau_m
	n' = (n_inf - n)/tau_n
	h' = (h_inf - h)/tau_h
}

FUNCTION xinf(v, v_half, k) { xinf = 1./(1+exp((v_half - v)/k)) }
FUNCTION tau(v, c_base, c_amp, v_max, sigma) { tau = c_base + c_amp*exp(-(v_max - v)*(v_max - v)/(sigma*sigma)) }

