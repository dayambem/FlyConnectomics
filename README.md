# Path Conductance

### Computational Connectomics Toolkit for Drosophila Neural Pathway Analysis

**PathConductance** is an R-based toolkit for identifying and ranking signal pathways in *Drosophila melanogaster* connectomic datasets. It integrates `neuprintr`, `natverse`, and custom graph algorithms to evaluate synaptic connectivity strength using a biologically-informed adaptation of Dijkstra’s algorithm.

---

## Features

- **Automated Neural Circuit Discovery** – Identify downstream neuron targets across multiple synaptic layers.
- **Path Conductance Metric** – Quantify connectivity using series and parallel synaptic resistance models.
- **Graph-Based Visualization** – Generate interactive hierarchical connectivity maps with `visNetwork`.
- **Heatmap-Based Connectivity Strength** – Visualize synaptic weight distributions with `pheatmap`.
- **Flexible Neuron Filtering** – Select neurons based on ROI, neurotransmitter type, or connection strength.

---

## How It Works

### Path Conductance Algorithm

Path conductance quantifies how efficiently signals propagate through neural circuits. Synaptic weights are treated as conductances (i.e., inverses of resistance). For a path A → B → D:

r_ABD = (1 / w_AB) + (1 / w_BD)

where `w_AB` and `w_BD` are synaptic weights. Local conductance is:

g_ABD_local = 1 / r_ABD

To compute global conductance over multiple parallel paths:

g_global = g_P1 + g_P2 + ... + g_Pn


This favors short, high-weight, and parallel routes—aligning with biological signal efficiency.

---

## Example Usage

```r
# Load dependencies and authenticate
library(neuprintr)
neuprint_login(server="https://neuprint.janelia.org", token="your_token")
source("path_conductance_functions.R")

# Get downstream neurons
downstream <- getDownstreamsNR(ids = c(123456), nLayers = 2, connStrength = 5)

# Compute path conductance matrix
cond_matrix <- computeConductance(downstream)

# Visualize network
plotNetworkGraph(cond_matrix)




