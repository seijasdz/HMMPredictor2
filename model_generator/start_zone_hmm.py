import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
from converter_to import converter_to
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import equal_distribution
from matrix_from_aln import matrix_from_exa


def test(model):
    with open('/home/zippyttech/Projects/personal/PerKalk/train_start2.exa') as in_file:
        oks = 0
        not_ok = 0
        for x_line in in_file:
            test_line = x_line.lower().replace('\n', '').replace(' ', '')
            tonight = converter_to(test_line, 2)
            logp, path = model.viterbi(tonight)
            path = [x[1].name for i, x in enumerate(path) if i < len(tonight)]
            if path[48] == 'start zone7':
                oks += 1
            else:
                not_ok += 1
        print(oks / (oks + not_ok))


back = State(DiscreteDistribution(equal_distribution), name='back')
back2 = State(DiscreteDistribution(equal_distribution), name='back2')

matrixZE = numpy.array(matrix_from_exa('new_tss.exa'))
start_states_data = classify(matrixZE, 2)
start_states = sequence_state_mulfactory(start_states_data, 'start zone')


model = HiddenMarkovModel()

model.add_state(back)
model.add_state(back2)
add_sequence(model, start_states)

model.add_transition(model.start, back, 1)
model.add_transition(back, back, 0.55)
model.add_transition(back, start_states[0], 0.45)
model.add_transition(start_states[-1], back2, 1)
model.add_transition(back2, back2, 0.5)


model.bake()

test(model)

lines = []
with open('/home/zippyttech/Projects/personal/PerKalk/train_start2.exa') as fi:
    for line in fi:
        lines.append(converter_to(line.replace('\n', '')))

model.fit(lines,
          transition_pseudocount=1,
          emission_pseudocount=1,
          verbose=True)

test(model)

with open('partial_model_start_model.json', 'w') as out:
    out.write(model.to_json())
