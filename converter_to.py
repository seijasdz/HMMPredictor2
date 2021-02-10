
def converter_to(sequence, order=2):
    new_list = []
    for index, element in enumerate(sequence):
        if index > order - 1:
            emission_name = element + '|'
            for x in range(order, 0, -1):
                emission_name += sequence[index - x]
            new_list.append(emission_name)
    return new_list


def converter_to2(seq):
    new = ''
    for i, e in enumerate(seq):
        if i > 1:
            emission = e + '|' + seq[i - 2] + seq[i - 1] + ','
            new += emission
    return new