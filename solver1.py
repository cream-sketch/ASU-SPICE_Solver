import numpy as np


def parse_spice_netlist(filename):

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
                # R<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('R', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)
            elif component_type == 'I':
                # I<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('I', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)
            elif component_type == 'V':
                # V<name> <node1> <node2> <value>
                node1, node2, value = parts[1], parts[2], float(parts[3])
                components.append(('V', node1, node2, value))
                nodes.add(node1)
                nodes.add(node2)


    if '0' in nodes:
        nodes.remove('0')
    nodes = sorted(list(nodes))
    return components, nodes


def build_system(components, nodes):
    """ G & J"""
    n = len(nodes)  # 节点数量
    num_voltage_sources = sum(1 for c in components if c[0] == 'V')
    total_unknowns = n + num_voltage_sources  # 总未知量（节点电压 + 电压源电流）

    G = np.zeros((total_unknowns, total_unknowns))
    J = np.zeros(total_unknowns)

    # 节点到索引的映射
    node_index = {node: idx for idx, node in enumerate(nodes)}

    # 电压源计数器
    vs_index = n  # 电压源电流的索引从 n 开始

    for component in components:
        c_type = component[0]

        if c_type == 'R':
            # 电阻
            _, node1, node2, value = component
            g = 1 / value  # 电导
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
            # 电流源
            _, node1, node2, value = component
            if node1 != '0':
                i = node_index[node1]
                J[i] -= value  # 流出
            if node2 != '0':
                j = node_index[node2]
                J[j] += value  # 流入

        elif c_type == 'V':
            # 电压源
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
    """求解电路的节点电压"""
    components, nodes = parse_spice_netlist(filename)
    G, J, nodes = build_system(components, nodes)

    # 求解 G * V = J
    V = np.linalg.solve(G, J)

    # 输出节点电压
    print(f"\nResults for {filename}:")
    for i, node in enumerate(nodes):
        print(f"Voltage at {node}: {V[i]:.3f} V")


# 测试两个电路
solve_circuit('circuit1.sp')
solve_circuit('circuit2.sp')
solve_circuit('circuit3.sp')