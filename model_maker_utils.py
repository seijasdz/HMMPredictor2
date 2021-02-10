from pomegranate import DiscreteDistribution
from pomegranate import State
from converter_to import converter_to

equal_distribution = {
    'a|aa': 1/64,
    'a|ac': 1/64,
    'a|ag': 1/64,
    'a|at': 1/64,
    'a|ca': 1/64,
    'a|cc': 1/64,
    'a|cg': 1/64,
    'a|ct': 1/64,
    'a|ga': 1/64,
    'a|gc': 1/64,
    'a|gg': 1/64,
    'a|gt': 1/64,
    'a|ta': 1/64,
    'a|tc': 1/64,
    'a|tg': 1/64,
    'a|tt': 1/64,
    'c|aa': 1/64,
    'c|ac': 1/64,
    'c|ag': 1/64,
    'c|at': 1/64,
    'c|ca': 1/64,
    'c|cc': 1/64,
    'c|cg': 1/64,
    'c|ct': 1/64,
    'c|ga': 1/64,
    'c|gc': 1/64,
    'c|gg': 1/64,
    'c|gt': 1/64,
    'c|ta': 1/64,
    'c|tc': 1/64,
    'c|tg': 1/64,
    'c|tt': 1/64,
    'g|aa': 1/64,
    'g|ac': 1/64,
    'g|ag': 1/64,
    'g|at': 1/64,
    'g|ca': 1/64,
    'g|cc': 1/64,
    'g|cg': 1/64,
    'g|ct': 1/64,
    'g|ga': 1/64,
    'g|gc': 1/64,
    'g|gg': 1/64,
    'g|gt': 1/64,
    'g|ta': 1/64,
    'g|tc': 1/64,
    'g|tg': 1/64,
    'g|tt': 1/64,
    't|aa': 1/64,
    't|ac': 1/64,
    't|ag': 1/64,
    't|at': 1/64,
    't|ca': 1/64,
    't|cc': 1/64,
    't|cg': 1/64,
    't|ct': 1/64,
    't|ga': 1/64,
    't|gc': 1/64,
    't|gg': 1/64,
    't|gt': 1/64,
    't|ta': 1/64,
    't|tc': 1/64,
    't|tg': 1/64,
    't|tt': 1/64
}

class HighOrderState:
    def __init__(self, columns):
        self.states_distribution = self.calculate_states(columns)

    @staticmethod
    def calculate_states(columns):

        states_distribution = {
            'a|aa': 1,
            'a|ac': 1,
            'a|ag': 1,
            'a|at': 1,
            'a|ca': 1,
            'a|cc': 1,
            'a|cg': 1,
            'a|ct': 1,
            'a|ga': 1,
            'a|gc': 1,
            'a|gg': 1,
            'a|gt': 1,
            'a|ta': 1,
            'a|tc': 1,
            'a|tg': 1,
            'a|tt': 1,
            'c|aa': 1,
            'c|ac': 1,
            'c|ag': 1,
            'c|at': 1,
            'c|ca': 1,
            'c|cc': 1,
            'c|cg': 1,
            'c|ct': 1,
            'c|ga': 1,
            'c|gc': 1,
            'c|gg': 1,
            'c|gt': 1,
            'c|ta': 1,
            'c|tc': 1,
            'c|tg': 1,
            'c|tt': 1,
            'g|aa': 1,
            'g|ac': 1,
            'g|ag': 1,
            'g|at': 1,
            'g|ca': 1,
            'g|cc': 1,
            'g|cg': 1,
            'g|ct': 1,
            'g|ga': 1,
            'g|gc': 1,
            'g|gg': 1,
            'g|gt': 1,
            'g|ta': 1,
            'g|tc': 1,
            'g|tg': 1,
            'g|tt': 1,
            't|aa': 1,
            't|ac': 1,
            't|ag': 1,
            't|at': 1,
            't|ca': 1,
            't|cc': 1,
            't|cg': 1,
            't|ct': 1,
            't|ga': 1,
            't|gc': 1,
            't|gg': 1,
            't|gt': 1,
            't|ta': 1,
            't|tc': 1,
            't|tg': 1,
            't|tt': 1
        }
        for index, base1 in enumerate(columns[0]):
            name = base1 + '|'
            for i, col in enumerate(columns):
                if i:
                    name += col[index]
            if name not in states_distribution:
                states_distribution[name] = 1
            else:
                states_distribution[name] += 1

        for key, value in states_distribution.items():
            states_distribution[key] /= len(columns[0])
        return states_distribution


def percentage_matrix_maker(seqs_data):
    matrix = []
    for data in seqs_data:
        for i in range(0, data[1]):
            matrix.append(list(data[0].lower()))
    return matrix


def sequence_state_factory(states_data, name):
    states = []
    for index, data in enumerate(states_data):
        state = State(DiscreteDistribution(data.states_distribution), name=name + str(index))
        states.append(state)
    return states


def spacer_states_maker(quantity, distribution, name):
    states = []
    for i in range(0, quantity):
        state = State(DiscreteDistribution(distribution), name=name + str(i))
        states.append(state)
    return states


def classify(matrix, order=1):
    state_data = []
    for index, column in enumerate(matrix.T):
        columns = [column]
        for x in range(order, 0, -1):
            col = matrix.T[index - x]
            columns.append(col)
        if index > order - 1:
            data = HighOrderState(columns)
            state_data.append(data)
    return state_data


def add_sequence(model, states):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 1.0)


def add_variable_length_sequence(model, states, end_state):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 0.75)
            model.add_transition(state, end_state, 0.25)
    model.add_transition(states[-1], end_state, 1.0)


def load_long_training_examples(filename, n):
    with open(filename) as file:
        first = file.read().replace('\n', '').replace(' ', '').split('>')[1:n]
        second = [converter_to(f.split('.')[1].lower(), 2) for f in first]
    return second


class StateNotFoundException(Exception):
    pass


def get_state(model, name):
    for state in model.states:
        if state.name == name:
            return state
    raise StateNotFoundException('State not found ' + name)


if __name__ == '__main__':
    seqs = [(['a', 'a'], 90), (['a', 'b'], 10)]
    print(percentage_matrix_maker(seqs))
