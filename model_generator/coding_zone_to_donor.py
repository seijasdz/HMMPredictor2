import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
import calculator
from converter_to import converter_to
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import equal_distribution
from matrix_from_aln import matrix_from_exa


def test_model(model):
    max = 0
    with open('new_ccutsa.txt') as ifile:
        examples = 0
        oks = 0
        for line in ifile:
            splitted = line.lower().split()[1:-1]
            for sp in splitted:
                examples += 1
                p_index = sp.find('p')

                test_line = converter_to(sp.replace('p', ''), 2)

                if len(test_line) > max:
                    max = len(test_line)

                logp, path = model.viterbi(test_line)
                path_names = [x[1].name for x in path]
                print('p', p_index)
                print('d', path_names.index('donor00'))
                for i, pn in enumerate(path_names):
                    if pn == 'donor00':
                        if p_index == i + 5:

                            oks +=1
        print('RESULTADO', oks / examples, 'max_len', max)

with open('new_cutsa.txt') as in_file:
    total = []
    for line in in_file:
        no_p_line = line.replace('P', '').lower().split()[:-1]
        total += no_p_line
converted_total = [converter_to(x, 2) for x in total]

matrixDonor0 = numpy.array(matrix_from_exa('new_donor1.exa'))

c0, c1, c2 = calculator.calculate_proba2('cuts.txt')

coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

donor0_data = classify(matrixDonor0, 2)
donor0_states = sequence_state_factory(donor0_data, 'donor0')

post = State(DiscreteDistribution(equal_distribution), name='post')

model = HiddenMarkovModel('codiing to donor')

model.add_state(coding_state0)
model.add_state(coding_state1)
model.add_state(coding_state2)

add_sequence(model, donor0_states)

model.add_state(post)

model.add_transition(model.start, coding_state0, 0.34)
model.add_transition(model.start, coding_state1, 0.33)
model.add_transition(model.start, coding_state2, 0.33)

model.add_transition(coding_state0, coding_state1,      0.6)
model.add_transition(coding_state0, donor0_states[0],   0.4)

model.add_transition(coding_state1, coding_state2,      0.6)
model.add_transition(coding_state1, donor0_states[0],   0.4)

model.add_transition(coding_state2, coding_state0,      0.6)
model.add_transition(coding_state2, donor0_states[0],   0.4)

model.add_transition(donor0_states[-1], post, 1)

model.add_transition(post, post, 0.9)
model.add_transition(post, model.end, 0.1)

model.bake()
test_model(model)

model.fit(converted_total,
          transition_pseudocount=1,
          emission_pseudocount=1,
          verbose=True)

test_model(model)


with open('partial_model_coding_to_donor_model.json', 'w') as out:
    out.write(model.to_json())