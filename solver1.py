import numpy as np

def parse_spice_netlist(filename):
    """Parse a SPICE netlist file and return the list of components and set of nodes."""
    components = []
    nodes = set()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*'):
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
    """Build the G matrix and J vector."""
    n = len(nodes)  # Number of nodes
    num_voltage_sources = sum(1 for c in components if c[0] == 'V')  # Count voltage sources
    total_unknowns = n + num_voltage_sources  # Total unknowns (node voltages + voltage source currents)

    G = np.zeros((total_unknowns, total_unknowns))  # Initialize G matrix
    J = np.zeros(total_unknowns)  # Initialize J vector

    # Map nodes to indices
    node_index = {node: idx for idx, node in enumerate(nodes)}

    # Counter for voltage source currents
    vs_index = n  # Index for voltage source currents starts at n

    for component in components:
        c_type = component[0]

        if c_type == 'R':
            # Resistor
            _, node1, node2, value = component
            g = 1 / value  # Conductance
            if node1 != '0' and node2 != '0':
                i, j = node_index[node1], node_index[node2]
                G[i][i] += g
                G[j][j] += g
                G[i][j] -= g
                G[j][i] -= g
            elif node1 != '0':
                i = node_index[node1]
                G[i][i] += g
            elif node2 != '0':
                j = node_index[node2]
                G[j][j] += g

        elif c_type == 'I':
            # Current source
            _, node1, node2, value = component
            if node1 != '0':
                i = node_index[node1]
                J[i] -= value  # Outflow
            if node2 != '0':
                j = node_index[node2]
                J[j] += value  # Inflow

        elif c_type == 'V':
            # Voltage source
            _, node1, node2, value = component
            if node1 != '0':
                i = node_index[node1]
                G[i][vs_index] = 1
                G[vs_index][i] = 1
            if node2 != '0':
                j = node_index[node2]
                G[j][vs_index] = -1
                G[vs_index][j] = -1
            J[vs_index] = value
            vs_index += 1

    return G, J, nodes

def solve_circuit(filename):
    """Solve for node voltages in the circuit."""
    components, nodes = parse_spice_netlist(filename)
    G, J, nodes = build_system(components, nodes)

    # Solve the system G * V = J
    V = np.linalg.solve(G, J)

    # Print node voltages
    print(f"\nResults for {filename}:")
    for i, node in enumerate(nodes):
        print(f"Voltage at {node}: {V[i]:.3f} V")

# Test the circuits
solve_circuit('circuit1.sp')
solve_circuit('circuit2.sp')