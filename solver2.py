import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg
import os

def parse_spice_netlist(filename):
    """Parse a SPICE netlist file and return the list of components and set of nodes."""
    components = []
    nodes = set()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*') or line.startswith('.'):
                continue
            parts = line.split()
            component_type = parts[0][0].upper()
            if component_type == 'R':
                # Resistor: R<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('R', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)
            elif component_type == 'I':
                # Current source: I<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('I', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)
            elif component_type == 'V':
                # Voltage source: V<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('V', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)

    # Remove the ground node '0'
    if '0' in nodes:
        nodes.remove('0')
    nodes = sorted(list(nodes))  # Sort nodes by name
    return components, nodes

def build_system(components, nodes):
    """Build the sparse G matrix and J vector."""
    # Identify nodes fixed by voltage sources
    fixed_nodes = {}
    for component in components:
        if component[0] == 'V':
            _, node1, node2, value = component
            if node2 == '0':  # Assume the negative terminal of the voltage source is connected to ground
                fixed_nodes[node1] = float(value)
            elif node1 == '0':
                fixed_nodes[node2] = -float(value)

    # Remove fixed nodes from active nodes
    active_nodes = [node for node in nodes if node not in fixed_nodes]
    n = len(active_nodes)
    # Correctly calculate the number of voltage sources that require additional unknowns
    num_voltage_sources = 0
    for c in components:
        if c[0] == 'V':
            _, node1, node2, _ = c
            # Only count voltage sources where at least one node is not fixed and not ground
            if (node1 not in fixed_nodes or node2 not in fixed_nodes) and node1 != '0' and node2 != '0':
                num_voltage_sources += 1
    total_unknowns = n + num_voltage_sources

    # Use sparse matrix (CSR format)
    row = []
    col = []
    data = []
    J = np.zeros(total_unknowns)

    node_index = {node: idx for idx, node in enumerate(active_nodes)}
    vs_index = n

    for component in components:
        c_type = component[0]

        if c_type == 'R':
            _, node1, node2, value = component
            g = 1 / value  # Conductance
            # Skip if both nodes are fixed or one node is ground and the other is fixed
            if (node1 in fixed_nodes and node2 in fixed_nodes) or (node1 == '0' and node2 in fixed_nodes) or (node2 == '0' and node1 in fixed_nodes):
                continue
            # node1 is fixed or ground, node2 is not fixed
            elif node1 in fixed_nodes or node1 == '0':
                if node2 not in fixed_nodes and node2 != '0':  # Ensure node2 is in active_nodes
                    j = node_index[node2]
                    row.append(j)
                    col.append(j)
                    data.append(g)
                    if node1 in fixed_nodes:
                        J[j] += g * fixed_nodes[node1]
                    # If node1 is '0', voltage is 0, no need to adjust J
            # node2 is fixed or ground, node1 is not fixed
            elif node2 in fixed_nodes or node2 == '0':
                if node1 not in fixed_nodes and node1 != '0':  # Ensure node1 is in active_nodes
                    i = node_index[node1]
                    row.append(i)
                    col.append(i)
                    data.append(g)
                    if node2 in fixed_nodes:
                        J[i] += g * fixed_nodes[node2]
                    # If node2 is '0', voltage is 0, no need to adjust J
            # Both nodes are not fixed and not ground
            else:
                i, j = node_index[node1], node_index[node2]
                # G[i][i] += g
                row.append(i)
                col.append(i)
                data.append(g)
                # G[j][j] += g
                row.append(j)
                col.append(j)
                data.append(g)
                # G[i][j] -= g
                row.append(i)
                col.append(j)
                data.append(-g)
                # G[j][i] -= g
                row.append(j)
                col.append(i)
                data.append(-g)

        elif c_type == 'I':
            _, node1, node2, value = component
            # Skip if both nodes are fixed or ground
            if (node1 in fixed_nodes and node2 in fixed_nodes) or (node1 == '0' and node2 in fixed_nodes) or (node2 == '0' and node1 in fixed_nodes) or (node1 == '0' and node2 == '0'):
                continue
            # node1 is fixed or ground, node2 is not fixed
            elif node1 in fixed_nodes or node1 == '0':
                if node2 not in fixed_nodes and node2 != '0':
                    j = node_index[node2]
                    J[j] += value
            # node2 is fixed or ground, node1 is not fixed
            elif node2 in fixed_nodes or node2 == '0':
                if node1 not in fixed_nodes and node1 != '0':
                    i = node_index[node1]
                    J[i] -= value
            # Both nodes are not fixed and not ground
            else:
                i, j = node_index[node1], node_index[node2]
                J[i] -= value
                J[j] += value

        elif c_type == 'V':
            _, node1, node2, value = component
            # Skip if both nodes are fixed or ground
            if (node1 in fixed_nodes or node1 == '0') and (node2 in fixed_nodes or node2 == '0'):
                continue
            # Only add an extra unknown if at least one node is not fixed and not ground
            if (node1 not in fixed_nodes or node2 not in fixed_nodes) and node1 != '0' and node2 != '0':
                if node1 not in fixed_nodes and node1 != '0':
                    i = node_index[node1]
                    row.append(i)
                    col.append(vs_index)
                    data.append(1)
                    row.append(vs_index)
                    col.append(i)
                    data.append(1)
                if node2 not in fixed_nodes and node2 != '0':
                    j = node_index[node2]
                    row.append(j)
                    col.append(vs_index)
                    data.append(-1)
                    row.append(vs_index)
                    col.append(j)
                    data.append(-1)
                J[vs_index] = value
                vs_index += 1

    # Build the sparse matrix
    G = sparse.csr_matrix((data, (row, col)), shape=(total_unknowns, total_unknowns))
    return G, J, active_nodes, fixed_nodes

def solve_circuit(filename):
    """Solve for node voltages in the circuit."""
    try:
        components, nodes = parse_spice_netlist(filename)
        G, J, active_nodes, fixed_nodes = build_system(components, nodes)

        print(f"\nNumber of non-zero elements in G matrix ({filename}): {G.nnz}")
        print(f"J vector ({filename}):")
        print(J)

        # Solve using sparse matrix solver
        V = splinalg.spsolve(G, J)

        # Combine fixed nodes and solved node voltages
        voltages = fixed_nodes.copy()
        for i, node in enumerate(active_nodes):
            voltages[node] = V[i]

        print(f"\nResults for {filename}:")
        for node in sorted(voltages.keys()):
            print(f"Voltage at {node}: {voltages[node]:.3f} V")

        return voltages
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None

def solve_all_circuits(directory, output_file):
    """Recursively traverse the directory, solve all .sp files, and save results to a file."""
    with open(output_file, 'w', encoding='utf-8') as f:
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if filename.endswith('.sp'):
                    filepath = os.path.join(root, filename)
                    f.write(f"\nResults for {filepath}:\n")
                    voltages = solve_circuit(filepath)
                    if voltages:
                        for node in sorted(voltages.keys()):
                            f.write(f"Voltage at {node}: {voltages[node]:.3f} V\n")

# Use full paths for the files
circuit1_path = r"G:\Root\PhD\ASU\assignment_simulation\circuit1.sp"
circuit2_path = r"G:\Root\PhD\ASU\assignment_simulation\circuit2.sp"

# Test the two circuits
solve_circuit(circuit1_path)
solve_circuit(circuit2_path)

# Traverse the GitHub dataset directory
directory = r"G:\Root\PhD\ASU\assignment_simulation\real-circuit-data"
output_file = r"G:\Root\PhD\ASU\assignment_simulation\results.txt"
solve_all_circuits(directory, output_file)