import galois
import math
import random
import copy


def initialize_field(k):
    field_order = 1 << math.ceil(math.log2(2 * math.ceil(math.log2(k)) + 2))
    field = galois.GF(field_order)
    return field, field_order


def find_primitive_root(GF):
    return GF.primitive_element


def e_poly(y, beta, field):
    res = field(1)
    for yi, betai in zip(y, beta):
        res *= field(1) - yi + (field(2) * yi - field(1)) * betai
    return res


def q_u_i(y, z, f_u, field, i):
    res = field(0)
    l = k.bit_length()
    for alpha in range(k):
        for beta in range(k):
            for gamma in range(k):
                alpha_bits = [int(bit) for bit in f'{alpha:0{l}b}'[::-1]]
                beta_bits = [int(bit) for bit in f'{beta:0{l}b}'[::-1]]
                gamma_bits = [int(bit) for bit in f'{gamma:0{l}b}'[::-1]]
                if alpha_bits[i] == 1 and f_u[beta][gamma] == alpha:
                    res += e_poly(y, beta_bits, field) * e_poly(z, gamma_bits, field)
    return res



def field_sum(bits, field):
    total = field(0)
    for i, bit in enumerate(bits):
        total += bit * (field(2) ** i)
    return total


def recursive_clean(node, level, h):
    log_k = math.ceil(math.log2(k))
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
            registers[node][i] -= q_u_i(registers[left_child], registers[right_child], tree[node], field, i)

        recursive_clean_inverse(left_child, level + 1, h)
        recursive_clean_inverse(right_child, level + 1, h)

        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] /= omega ** j


def recursive_clean_inverse(node, level, h):
    log_k = math.ceil(math.log2(k))
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

        recursive_clean(left_child, level + 1, h)
        recursive_clean(right_child, level + 1, h)

        for i in range(log_k):
            registers[node][i] += q_u_i(registers[left_child], registers[right_child], tree[node], field, i)

        recursive_clean_inverse(left_child, level + 1, h)
        recursive_clean_inverse(right_child, level + 1, h)

        for c in [left_child, right_child]:
            for i in range(log_k):
                registers[c][i] /= omega ** j



def clean_computation(h):
    initial_registers = copy.deepcopy(registers[0])
    recursive_clean(0, 0, h)
    final_result = field_sum(registers[0], field) - field_sum(initial_registers, field)
    return final_result, registers


def initialize_tree_and_catalyst():
    h = 1  # int(input("Enter tree height: "))
    k = 10  # int(input("Enter alphabet size: "))
    num_nodes = (1 << (h + 1)) - 1

    tree_values = [1] * (1 << h) # to be defined
    tree_functions = [[list(range(k)) for _ in range(k)] for _ in range((1 << h) - 1)]  # to be defined
    tree = tree_functions + tree_values

    field, m = initialize_field(k)

    registers = [[field(random.randint(0, field.order - 1)) for _ in range(math.ceil(math.log2(k)))] for _ in
                range(num_nodes)]

    return h, k, tree, registers, m, field


if __name__ == "__main__":
    height, k, tree, registers, m, field = initialize_tree_and_catalyst()
    omega = find_primitive_root(field)
    registers_copy = copy.deepcopy(registers)
    result, final_registers = clean_computation(height)
    print(f"Tree evaluation result: {result}")
    print(f"Initial registers: {registers_copy}")
    print(f"Final registers: {final_registers}")
