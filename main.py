import galois
import math
import random
import copy


def q_u_i(y, z, f_u, i):
    def e_poly(y, beta):
        res = field(1)
        for yi, betai in zip(y, beta):
            res *= field(1) - yi + (2 * yi - field(1)) * betai
        return res

    res = field(0)
    for alpha in range(k):
        for beta in range(k):
            for gamma in range(k):
                alpha_bits = [field(int(bit)) for bit in f'{alpha:0{log_k}b}'[::-1]]
                beta_bits = [field(int(bit)) for bit in f'{beta:0{log_k}b}'[::-1]]
                gamma_bits = [field(int(bit)) for bit in f'{gamma:0{log_k}b}'[::-1]]
                if alpha_bits[i] == field(1) and f_u[beta][gamma] == alpha:
                    res += e_poly(y, beta_bits) * e_poly(z, gamma_bits)
    return res


def field_sum(bits):
    total = field(0)
    for i, bit in enumerate(bits):
        total += bit * (field(2) ** i)
    return total


def recursive_clean(node, level, h):
    if level == h:
        bits = [field(int(bit)) for bit in f'{tree[node]:0{log_k}b}'[::-1]]
        for i, bit in enumerate(bits):
            registers[node][i] += bit
        return

    left_child = 2 * node + 1
    right_child = 2 * node + 2

    for j in range(1, m + 1):
        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] *= omega ** j

        recursive_clean(left_child, level + 1, h)
        recursive_clean(right_child, level + 1, h)

        for i in range(log_k):
            registers[node][i] -= q_u_i(registers[left_child], registers[right_child], tree[node], i)

        recursive_clean_inverse(left_child, level + 1, h)
        recursive_clean_inverse(right_child, level + 1, h)

        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] /= omega ** j


def recursive_clean_inverse(node, level, h):
    if level == h:
        bits = [field(int(bit)) for bit in f'{tree[node]:0{log_k}b}'[::-1]]
        for i, bit in enumerate(bits):
            registers[node][i] -= bit
        return

    left_child = 2 * node + 1
    right_child = 2 * node + 2

    for j in range(1, m + 1):
        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] *= omega ** j

        recursive_clean_inverse(left_child, level + 1, h)
        recursive_clean_inverse(right_child, level + 1, h)

        for i in range(log_k):
            registers[node][i] += q_u_i(registers[left_child], registers[right_child], tree[node], i)

        recursive_clean(left_child, level + 1, h)
        recursive_clean(right_child, level + 1, h)

        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] /= omega ** j


def clean_computation(h):
    initial_registers = copy.deepcopy(registers[0])
    recursive_clean(0, 0, h)
    final_result = field_sum(registers[0]) - field_sum(initial_registers)
    return final_result, registers


def initialize_field(k):
    field_order = 1 << math.ceil(math.log2(2 * math.ceil(math.log2(k)) + 2))
    field = galois.GF(field_order)
    return field, field_order


def initialize_tree_and_catalyst():
    h = int(input("Enter the height of the tree: "))
    k = int(input("Enter the size of the alphabet: "))
    num_nodes = (1 << (h + 1)) - 1
    function = [[int(input(f"Enter the value for function at row {r}, column {c}: ")) for c in range(k)]
                for r in range(k)]
    tree_functions = [function for _ in range((1 << h) - 1)]

    tree_values = [int(input(f"Enter the value of leaf node {i}: ")) for i in range(1 << h)]
    tree = tree_functions + tree_values

    log_k = math.ceil(math.log2(k + 1))

    field, field_order = initialize_field(k)
    registers = [[field(random.randint(0, field.order - 1)) for _ in range(log_k)] for _ in
                 range(num_nodes)]

    return h, k, tree, registers, field, field_order - 1


def the_not_so_cool_algorithm(node, level):
    if level == height:
        return tree[node]

    left = the_not_so_cool_algorithm(2 * node + 1, level + 1)
    right = the_not_so_cool_algorithm(2 * node + 2, level + 1)

    return tree[node][left][right]


if __name__ == "__main__":
    height, k, tree, registers, field, m = initialize_tree_and_catalyst()
    log_k = math.ceil(math.log2(k + 1))
    omega = field.primitive_element
    registers_copy = copy.deepcopy(registers)
    result, final_registers = clean_computation(height)

    print("*" * 80)
    assert result == the_not_so_cool_algorithm(0, 0)
    print(f"Tree evaluation result: {result}")
    print(f"Initial registers:\n{registers_copy}")
    print(f"Final registers:\n{final_registers}")
