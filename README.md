# 🧪 Advanced Chemical Equation Balancer & Simulator

## 📌 Project Overview
This project balances chemical equations using concepts from linear algebra. It converts chemical reactions into a system of linear equations and solves them using matrix operations such as RREF and null space.

An interactive interface is built using Streamlit to visualize each step of the solution.

---

## 🚀 Features

- ✅ Automatic chemical equation balancing  
- ✅ Matrix representation of equations  
- ✅ Step-by-step solution (RREF, rank, null space)  
- ✅ Detection of invalid/impossible equations  
- ✅ Mole ratio visualization (graph)  
- ✅ Limiting reagent simulator  

---

## 🧠 Concepts Used

- Systems of Linear Equations  
- Matrix Representation  
- Row Reduced Echelon Form (RREF)  
- Rank and Nullity  
- Null Space  

---

## 💻 Tech Stack

- Python  
- Streamlit  
- SymPy  
- Pandas  
- Matplotlib  

---

## ▶️ How to Run

1. Install dependencies:
```bash
pip install streamlit sympy pandas matplotlib
```
2. Run the app:
streamlit run app.py

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

## 📊 Project Workflow

1. Input chemical equation  
2. Convert into matrix form  
3. Solve Ax = 0  
4. Apply RREF  
5. Find null space  
6. Scale to integers  
7. Display balanced equation  

---

## 🎯 Applications

- Chemistry education  
- Demonstration of linear algebra concepts  
- Interactive learning tool  

---

## 👨‍💻 Author

- Pallavi B