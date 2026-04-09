import streamlit as st
import re
import pandas as pd
import matplotlib.pyplot as plt
from sympy import Matrix, lcm

st.title("🧪 Advanced Chemical Equation Balancer & Simulator")

# -------------------------------
# Parse compound
# -------------------------------
def parse_compound(compound):
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', compound)
    comp_dict = {}
    for el, count in elements:
        count = int(count) if count else 1
        comp_dict[el] = comp_dict.get(el, 0) + count
    return comp_dict


# -------------------------------
# Build matrix
# -------------------------------
def build_matrix(reactants, products):
    compounds = reactants + products
    elements = set()

    parsed_all = []
    for c in compounds:
        parsed = parse_compound(c)
        parsed_all.append(parsed)
        elements.update(parsed.keys())

    elements = list(elements)
    matrix = []

    for el in elements:
        row = []
        for c in reactants:
            row.append(parse_compound(c).get(el, 0))
        for c in products:
            row.append(-parse_compound(c).get(el, 0))
        matrix.append(row)

    return Matrix(matrix), elements


# -------------------------------
# Check element mismatch
# -------------------------------
def check_elements(reactants, products):
    left_elements = set()
    right_elements = set()

    for c in reactants:
        left_elements.update(parse_compound(c).keys())

    for c in products:
        right_elements.update(parse_compound(c).keys())

    return left_elements == right_elements


# -------------------------------
# Solve using null space
# -------------------------------
def solve(matrix):
    ns = matrix.nullspace()
    sol = ns[0]

    lcm_val = lcm([term.q for term in sol])
    coeffs = [int(term * lcm_val) for term in sol]

    return coeffs, sol, ns


# -------------------------------
# Format equation
# -------------------------------
def format_eq(reactants, products, coeffs):
    r_len = len(reactants)
    left = " + ".join(f"{coeffs[i]}{reactants[i]}" for i in range(r_len))
    right = " + ".join(f"{coeffs[i+r_len]}{products[i]}" for i in range(len(products)))
    return left + " -> " + right


# -------------------------------
# UI Input
# -------------------------------
equation = st.text_input("Enter equation (e.g., C2H6 + O2 -> CO2 + H2O)")

if "run" not in st.session_state:
    st.session_state.run = False

if st.button("Balance & Analyze"):
    st.session_state.run = True

if st.session_state.run:

    try:
        left, right = equation.split("->")
        reactants = [x.strip() for x in left.split("+")]
        products = [x.strip() for x in right.split("+")]

        # -------------------------------
        # Edge Case: Element mismatch
        # -------------------------------
        if not check_elements(reactants, products):
            st.error("❌ Invalid equation: Elements mismatch between reactants and products.")
            st.stop()

        # -------------------------------
        # Parsed compounds
        # -------------------------------
        st.subheader("🔷 Parsed Compounds")
        for comp in reactants + products:
            st.write(f"{comp} → {parse_compound(comp)}")

        # -------------------------------
        # Matrix
        # -------------------------------
        matrix, elements = build_matrix(reactants, products)

        st.subheader("🔷 Elements Identified")
        st.write(elements)

        df = pd.DataFrame(matrix.tolist(),index=elements,columns=reactants + products)

        st.subheader("🔷 Matrix Representation")
        st.dataframe(df)

        # -------------------------------
        # Equations
        # -------------------------------
        st.subheader("🔷 System of Equations")
        for i, el in enumerate(elements):
            eq = ""
            for j, comp in enumerate(reactants + products):
                coeff = matrix[i, j]
                if coeff != 0:
                    eq += f"{coeff}x{j+1} "
            eq += "= 0"
            st.write(f"{el}: {eq}")

        # -------------------------------
        # RREF
        # -------------------------------
        st.subheader("🔷 RREF")
        rref_matrix = matrix.rref()[0]
        st.write(rref_matrix)
        st.info("RREF removes redundant equations and shows independent constraints.")

        # -------------------------------
        # Rank & Nullity
        # -------------------------------
        st.subheader("🔷 Rank & Nullity")
        rank = matrix.rank()
        nullity = len(reactants + products) - rank
        st.write(f"Rank: {rank}")
        st.write(f"Nullity: {nullity}")

        # -------------------------------
        # Null space
        # -------------------------------
        st.subheader("🔷 Null Space")
        ns = matrix.nullspace()
        st.write(ns)

        # -------------------------------
        # Solve
        # -------------------------------
        coeffs, raw_sol, ns = solve(matrix)

        st.subheader("🔷 Raw Solution")
        st.write(raw_sol)

        st.subheader("🔷 Scaled Integer Solution")
        st.write(coeffs)

        # -------------------------------
        # Final equation
        # -------------------------------
        balanced = format_eq(reactants, products, coeffs)

        st.subheader("✅ Balanced Equation")
        st.success(balanced)

        # -------------------------------
        # Ratio Visualization
        # -------------------------------
        st.subheader("📊 Mole Ratio Visualization")

        labels = reactants + products
        values = coeffs

        fig, ax = plt.subplots()
        ax.bar(labels, values)
        ax.set_ylabel("Moles Ratio")
        st.pyplot(fig)

        # -------------------------------
        # Limiting Reagent
        # -------------------------------
        st.subheader("🧪 Limiting Reagent Simulator")

        user_inputs = []
        for i, r in enumerate(reactants):
            val = st.number_input(f"Moles of {r}", min_value=0.0, value=1.0)
            user_inputs.append(val)

        ratios = [user_inputs[i] / coeffs[i] for i in range(len(reactants))]
        limiting_index = ratios.index(min(ratios))

        st.write(f"Limiting Reagent: {reactants[limiting_index]}")

        # product formation
        for i in range(len(products)):
            produced = ratios[limiting_index] * coeffs[len(reactants) + i]
            st.write(f"{products[i]} produced: {produced:.2f} moles")

    except:
        st.error("Invalid input format. Use: A + B -> C")