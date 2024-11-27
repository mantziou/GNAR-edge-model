# GNAR-edge-model
Code and synthetic data experiments for the paper "The GNAR-edge model: A network autoregressive model for networks with time varying edge weights", https://doi.org/10.1093/comnet/cnad039

In script ar_ts_edges.R is the main code for the GNAR-edge model. Main code in script ar_ts_edges.R is modified version of GNAR model code by Knight et al. 2019, to address the problem when time series are observed on the edges of a network. In addition, code in ar_ts_edges.R includes option for lead-lag weights.

Script sim_gnar_edge_1st_sec.R contains code for simulation experiments of section 5.1 of the main paper, and script sim_gnar_edge_2nd_sec.R contains code for simulation experiments of section 5.2 of the main paper.
