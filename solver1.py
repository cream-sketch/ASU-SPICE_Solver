import numpy as np

def parse_spice_netlist(filename):
    """Parse a SPICE netlist file."""
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
    """Build the G matrix and J vector."""
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
    # Count voltage sources that require additional unknowns
    num_voltage_sources = 0
    for c in components:
        if c[0] == 'V':
            _, node1, node2, _ = c
            # Only count voltage sources where at least one node is not fixed and not ground
            if (node1 not in fixed_nodes and node2 not in fixed_nodes) and node1 != '0' and node2 != '0':
                num_voltage_sources += 1
    total_unknowns = n + num_voltage_sources

    # Initialize G matrix and J vector
    G = np.zeros((total_unknowns, total_unknowns))
    J = np.zeros(total_unknowns)

    # Map active nodes to indices
    node_index = {node: idx for idx, node in enumerate(active_nodes)}
    vs_index = n  # Index for voltage source currents starts at n

    for component in components:
        c_type = component[0]

        if c_type == 'R':
            # Resistor
            _, node1, node2, value = component
            g = 1 / value  # Conductance
            # Skip if both nodes are fixed or one node is ground and the other is fixed
            if (node1 in fixed_nodes and node2 in fixed_nodes) or (node1 == '0' and node2 in fixed_nodes) or (node2 == '0' and node1 in fixed_nodes):
                continue
            # node1 is fixed or ground, node2 is not fixed
            elif node1 in fixed_nodes or node1 == '0':
                if node2 not in fixed_nodes and node2 != '0':  # Ensure node2 is in active_nodes
                    j = node_index[node2]
                    G[j][j] += g
                    if node1 in fixed_nodes:
                        J[j] += g * fixed_nodes[node1]
                    # If node1 is '0', voltage is 0, no need to adjust J
            # node2 is fixed or ground, node1 is not fixed
            elif node2 in fixed_nodes or node2 == '0':
                if node1 not in fixed_nodes and node1 != '0':  # Ensure node1 is in active_nodes
                    i = node_index[node1]
                    G[i][i] += g
                    if node2 in fixed_nodes:
                        J[i] += g * fixed_nodes[node2]
                    # If node2 is '0', voltage is 0, no need to adjust J
            # Both nodes are not fixed and not ground
            else:
                i, j = node_index[node1], node_index[node2]
                G[i][i] += g
                G[j][j] += g
                G[i][j] -= g
                G[j][i] -= g

        elif c_type == 'I':
            # Current source
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
            # Voltage source
            _, node1, node2, value = component
            # Skip if both nodes are fixed or ground
            if (node1 in fixed_nodes or node1 == '0') and (node2 in fixed_nodes or node2 == '0'):
                continue
            # Only add an extra unknown if both nodes are not fixed and not ground
            if (node1 not in fixed_nodes and node2 not in fixed_nodes) and node1 != '0' and node2 != '0':
                i = node_index[node1]
                j = node_index[node2]
                G[i][vs_index] = 1
                G[vs_index][i] = 1
                G[j][vs_index] = -1
                G[vs_index][j] = -1
                J[vs_index] = value
                vs_index += 1

    return G, J, active_nodes, fixed_nodes

def solve_circuit(filename):
    """Solve for node voltages in the circuit."""
    components, nodes = parse_spice_netlist(filename)
    G, J, active_nodes, fixed_nodes = build_system(components, nodes)

    # Solve the system G * V = J
    V = np.linalg.solve(G, J)

    # Combine fixed nodes and solved node voltages
    voltages = fixed_nodes.copy()
    for i, node in enumerate(active_nodes):
        voltages[node] = V[i]

    # Print node voltages
    print(f"\nResults for {filename}:")
    for node in sorted(voltages.keys()):
        print(f"Voltage at {node}: {voltages[node]:.3f} V")

# Test the circuits
solve_circuit('circuit1.sp')
solve_circuit('circuit2.sp')