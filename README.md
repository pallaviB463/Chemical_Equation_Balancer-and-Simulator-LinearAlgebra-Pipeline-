# 🧪 Advanced Chemical Equation Balancer & Simulator

## 📌 Project Overview

This project balances chemical equations using concepts from linear algebra and network theory. It converts chemical reactions into a system of linear equations, solves them using matrix operations (RREF and null space), and analyzes thermodynamic feasibility.

Additionally, it builds a reaction network graph to find optimal pathways between compounds, detect thermodynamic cycles, and analyze network topology using graph algorithms (Dijkstra's, cycle detection, centrality measures).

An interactive interface built with Streamlit provides visualization for both equation balancing and reaction network analysis.

---

## 🚀 Features

- ✅ Automatic chemical equation balancing  
- ✅ Matrix representation of equations  
- ✅ Step-by-step solution (RREF, rank, null space)  
- ✅ Detection of invalid/impossible equations  
- ✅ Thermodynamic analysis (ΔG°f values)
- ✅ Reaction network graph visualization
- ✅ Shortest path finding between compounds (Dijkstra's algorithm)
- ✅ Cycle and loop detection in reaction networks
- ✅ Graph metrics (centrality, betweenness)
- ✅ Support for complex compounds with parentheses (e.g., Ca(OH)₂, Al₂(SO₄)₃)

---

## 🧠 Concepts Used

### Linear Algebra
- Systems of Linear Equations  
- Matrix Representation  
- Row Reduced Echelon Form (RREF)  
- Rank and Nullity  
- Null Space  

### Thermodynamics
- Gibbs Free Energy (ΔG°f)
- Hess's Law
- Reaction spontaneity

### Graph Theory & Algorithms
- Directed Weighted Multigraphs
- Dijkstra's Algorithm (shortest path)
- Cycle Detection (DFS)
- Tarjan's Algorithm (strongly connected components)
- Network centrality measures  

---

## 💻 Tech Stack

- Python  
- Streamlit  
- SymPy  
- NetworkX  
- Pandas  
- Matplotlib  
- NumPy  

---

## ▶️ How to Run

1. Install dependencies:
```bash
pip install streamlit sympy pandas matplotlib networkx numpy
```
2. Run the app:
```bash
streamlit run app.py
```

---

## 📸 Example

Input:
```bash
C2H6 + O2 -> CO2 + H2O
```
Output:
```bash
2C2H6 + 7O2 -> 4CO2 + 6H2O
```

---

## 🎯 Workflow

### Equation Balancer Tab
1. Input chemical equation  
2. Convert into matrix form  
3. Solve Ax = 0  
4. Apply RREF  
5. Find null space  
6. Scale to integers  
7. Display balanced equation  
8. Compute thermodynamic properties (ΔG°rxn)

### Reaction Network Graph Tab
1. Build directed weighted multigraph from reaction database
2. Visualize compound nodes and reaction edges
3. Find shortest thermodynamic path between compounds
4. Detect cycles and equilibrium loops
5. Analyze network topology and centrality metrics
6. Identify key intermediates and hub compounds  

---

## 🎯 Applications

- Chemistry education  
- Demonstration of linear algebra concepts  
- Network analysis of complex reaction pathways
- Thermodynamic feasibility assessment
- Reaction pathway optimization
- Interactive learning tool for chemistry and graph theory  

---

## 👨‍💻 Author

- Pallavi B
