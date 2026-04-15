# 🧪 Advanced Chemical Equation Balancer & Simulator

## 📌 Project Overview

A Streamlit web application that balances any chemical equation using **linear algebra** — specifically the null space of the element-conservation matrix. 

The app converts chemical reactions into a system of linear equations, solves them using matrix operations (RREF and null space), and presents a complete mathematical breakdown. It includes an interactive limiting reagent simulator to calculate consumption and product formation.

---

## 🚀 Features

| Feature | What it does | Why it matters |
|---|---|---|
| **Recursive formula parser** | Handles `Ca(OH)₂`, `Al₂(SO₄)₃`, `Fe₂(SO₄)₃` — any nested parentheses | Correctly expands complex compound formulas |
| **GCD reduction** | Scales by LCM then divides by GCD | Guarantees minimal integer coefficients — `[2,7,4,6]` not `[4,14,8,12]` |
| **Multi-solution handler** | Detects nullity > 1, enumerates all basis vectors | Handles equations with multiple independent solutions |
| **Precise mismatch reporting** | Shows left-only and right-only element sets | Tells you *exactly* which elements are missing on which side |

### Main Features
- ✅ Automatic chemical equation balancing  
- ✅ Element conservation matrix display
- ✅ RREF (Reduced Row Echelon Form)  
- ✅ Rank & Nullity analysis
- ✅ Null space computation (SymPy exact arithmetic)  
- ✅ Step-by-step solution breakdown  
- ✅ Detection of invalid/impossible equations  
- ✅ Balanced equation verification
- ✅ Mole ratio visualization (interactive chart)
- ✅ Limiting reagent simulator with excess calculation

---

## 🧠 Concepts Used

### Linear Algebra
- Systems of Linear Equations (Ax = 0)
- Matrix Representation  
- Row Reduced Echelon Form (RREF)  
- Rank and Nullity  
- Null Space  
- Exact rational arithmetic (SymPy vs NumPy)

### Stoichiometry & Chemistry
- Law of conservation of mass  
- Element balancing
- Coefficient ratios
- Limiting reagents
- Mole calculations  

---

## 📐 How It Works — The Mathematics

### Step 1: Parse the formula
The recursive parser reads each formula character by character. When it encounters `(`, it recurses into a sub-call and multiplies the returned counts by the subscript after `)`.

```
Ca(OH)2  →  Ca=1, then recurse into (OH): O=1, H=1, multiply by 2
         →  {Ca:1, O:2, H:2}
```

### Step 2: Build the matrix A
Each row = one element. Each column = one compound. Reactants → positive, products → negative. This encodes the law of conservation of mass as the homogeneous linear system **Ax = 0**.

```
C₂H₆ + O₂ → CO₂ + H₂O

     C₂H₆  O₂  CO₂  H₂O
C  [  2     0   -1    0  ]
H  [  6     0    0   -2  ]
O  [  0     2   -2   -1  ]
```

### Step 3: RREF and null space
SymPy computes the Reduced Row Echelon Form (exact rational arithmetic — no floating point errors). The null space basis vectors are read directly from RREF by setting free variables to 1.

```
RREF(A) = [[1, 0, 0, -1/3],
           [0, 1, 0, -7/6],
           [0, 0, 1, -2/3]]

Rank = 3, Nullity = 1
Null vector: [1/3, 7/6, 2/3, 1]
```

### Step 4: LCM scaling + GCD reduction
```
Denominators: [3, 6, 3, 1]   →   LCM = 6
Scale:  [1/3, 7/6, 2/3, 1] × 6 = [2, 7, 4, 6]
GCD(2, 7, 4, 6) = 1           →   already minimal
Result: 2C₂H₆ + 7O₂ → 4CO₂ + 6H₂O  ✓
```

### Why SymPy instead of NumPy?
SymPy uses **exact rational arithmetic**. `1/3` is stored as the fraction `Rational(1,3)`, not the float `0.333...`. This matters because `0.333... × 6 = 1.9999...` in floating point, which breaks integer conversion. SymPy guarantees `Rational(1,3) × 6 = 2` exactly.

---

## 💻 Tech Stack

- Python  
- Streamlit  
- SymPy (for exact rational arithmetic matrix operations)  
- Pandas (for data table display)
- Matplotlib (for charts)
- NumPy    

---

## 📦 Installation

```bash
pip install streamlit sympy pandas matplotlib numpy
```

---

## ▶️ How to Run

Run the app:
```bash
streamlit run app.py
```

---

## 💡 Usage

1. Type any equation in the input box, e.g.:
   - `C2H6 + O2 -> CO2 + H2O`
   - `Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O`
   - `Al2(SO4)3 + NaOH -> Al(OH)3 + Na2SO4`
2. Click **Balance ⚗️** or use a quick-example button
3. Scroll through the full mathematical breakdown

---

## 📋 Supported Formula Syntax

| Syntax | Example | Supported |
|---|---|---|
| Simple compounds | `H2O`, `CO2` | ✅ |
| Multi-letter elements | `Fe`, `Ca`, `Mg` | ✅ |
| Single parentheses | `Ca(OH)2` | ✅ |
| Nested parentheses | `Al2(SO4)3`, `Fe2(SO4)3` | ✅ |
| Decimal coefficients | `3.5 H2O` | ❌ (integers only) |

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

## 🎯 App Workflow

### Main Sections (Top to Bottom)

1. **Input** → Enter equation (e.g., `Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O`)
2. **Parsed Compounds** → Shows how each compound's formula is expanded  
3. **Element Conservation Matrix** → Rows = elements, Columns = compounds  
4. **System of Linear Equations** → Human-readable form (Ax = 0)
5. **RREF** → Reduced Row Echelon Form with pivot columns identified
6. **Rank & Nullity** → Determines number of independent solutions
7. **Null Space** → Basis vectors encoding all solutions
8. **Multi-Solution Selector** (if nullity > 1) → Choose which basis vector to use  
9. **Coefficient Calculation** → LCM scaling then GCD reduction to get minimal integers
10. **Balanced Equation** → Final verified result  
11. **Mole Ratio Chart** → Visual bar chart of stoichiometric coefficients  
12. **Limiting Reagent Simulator** → Enter reactant amounts, see limiting reagent and products formed

---

## 📚 Example Workflow

**Input:** `C2H6 + O2 -> CO2 + H2O`

**Output:** `2C2H6 + 7O2 → 4CO2 + 6H2O`

The app shows:
- How the recursive parser expands each compound
- The 3×4 conservation matrix for C, H, O
- The RREF computation  
- Why nullity = 1 (one unique solution)
- The LCM/GCD reduction steps
- A bar chart showing coefficients [2, 7, 4, 6]
- Interactive fields to input mole amounts and find which reactant runs out first

---

## 🎓 Educational Value

- Demonstrates core linear algebra: matrix rank, nullity, null space
- Shows practical application of exact arithmetic (SymPy vs NumPy)
- Teaches chemical stoichiometry interactively  
- Illustrates the recursive descent parsing algorithm
- Interactive limiting reagent analysis

---

## 👨‍💻 Author

- Pallavi B
