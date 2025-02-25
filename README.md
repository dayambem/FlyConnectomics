# FlyConnectomics
### **FlyConnectomics: Mapping Drosophila Neural Circuits with Network Science** 🧠🔬  

**FlyConnectomics** is an R-based computational toolkit designed to **analyze, predict, and visualize neural connectivity in Drosophila**. By leveraging **neuprintR, natverse, and graph theory**, this repository enables researchers to explore how **neural circuits drive complex behaviors** such as motor control, aggression, and flight.  

---

## **📌 Features**  

✅ **Automated Neural Circuit Discovery** – Identify downstream neuron targets across multiple synaptic layers.  
✅ **Path Inefficiency Analysis** – Quantify neural connectivity using a custom network science approach.  
✅ **Graph-Based Visualization** – Generate interactive **hierarchical connectivity maps** with **visNetwork**.  
✅ **Heatmap-Based Connectivity Strength** – Apply **pheatmap** to analyze synaptic weight distributions.  
✅ **Flexible Neuron Filtering** – Select neurons based on ROI, neurotransmitter type, or connectivity thresholds.  

---

## **📖 How It Works**  

### **1️⃣ Path Inefficiency Algorithm**  
Path inefficiency quantifies how **efficiently signals travel between neurons**. Path inefficiency is defined as:

r_ABD = r_AB + r_BD = (1 / w_AB) + (1 / w_BD)

where w_AB and w_BD are synaptic weights. The global inefficiency metric integrates all possible pathways for robust circuit mapping.  

### **2️⃣ Neural Circuit Identification**  
- **Extracts connectivity data** from [neuprint.janelia.org](https://neuprint.janelia.org).  
- **Applies shortest path algorithms** to predict functional neuron connections.  
- **Identifies strong downstream targets** of motor and descending neurons (DNs).  

### **3️⃣ Visualization & Analysis**  
- **Hierarchical graphs** depict neuron-to-motor pathways.  
- **Heatmaps** display synaptic strengths and connectivity inefficiencies.  
- **Graph-theoretic clustering** reveals key functional modules.  

---

## **🚀 Installation**  
Ensure you have **R (4.0+), RStudio, and the following packages installed:**  

```r
install.packages(c("neuprintr", "natverse", "igraph", "visNetwork", "pheatmap", "dplyr", "R.matlab"))
```
Clone the repository:  
```sh
git clone https://github.com/yourusername/FlyConnectomics.git
```

---

## **📊 Example Usage**  

### **Load Dependencies & Authenticate with NeuPrint**  
```r
library(natverse)
library(neuprintr)
neuprint_login(server="https://neuprint.janelia.org", token="your_token")
```

### **Retrieve & Analyze Neuronal Connections**  
```r
source("connectomics_functions.R")

# Get downstream neurons for a given input
downstream_neurons <- getDownstreamsNR(ids=c(123456), nLayers=2, connStrength=5)

# Compute path inefficiency for connections
efficiency_matrix <- computeInefficiency(downstream_neurons)
```

### **Visualize Neural Networks**  
```r
plotNetworkGraph(efficiency_matrix)
```

---

## **📎 Example Outputs**  

### **Global Path Inefficiency Heatmap**  
![Heatmap](image.png)  

### **Neural Circuit Visualization**  
![Network Graph](network_graph.png)  

---

## **🛠 Applications**  
✔ **Motor control & descending neuron analysis**  
✔ **Behavioral neuroscience (aggression, flight, learning circuits)**  
✔ **Graph-based connectomics modeling**  
✔ **Comparison with AI-driven neuroscience models**  

---

## **🤝 Contributing**  
We welcome **collaborations, issues, and pull requests**! Feel free to fork the repo and enhance the toolkit.  

---

## **📧 Contact**  
For questions, reach out via **[your email]** or visit [your website/linkedin].  

🚀 *Unlock the power of neural networks—one synapse at a time!*
