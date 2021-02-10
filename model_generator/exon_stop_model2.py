import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
import calculator
from stop_example_divider import divider as stop_divider
from model_maker_utils import sequence_state_factory, classify, add_sequence, equal_distribution
from matrix_from_aln import matrix_from_exa
from converter_to import converter_to

c0, c1, c2 = calculator.calculate_proba2('cuts.txt')
matrixStop = numpy.array(matrix_from_exa('new_tts.exa'))
coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

taa_matrix, tga_matrix, tag_matrix = stop_divider('new_tts.exa')

post = State(DiscreteDistribution(equal_distribution), name='post')

model = HiddenMarkovModel('coding_to_stop')

stop_data = classify(matrixStop, 2)
stop_states = sequence_state_factory(stop_data, 'stop')

model.add_state(coding_state0)
model.add_state(coding_state1)
model.add_state(coding_state2)

add_sequence(model, stop_states)

model.add_state(post)

model.add_transition(model.start, coding_state2, 1)
model.add_transition(coding_state0, coding_state1, 1)
model.add_transition(coding_state1, coding_state2, 1)
model.add_transition(coding_state2, coding_state0,  0.6)
model.add_transition(coding_state2, stop_states[0], 0.4)
model.add_transition(stop_states[-1], post, 1)
model.add_transition(post, post, 0.9)
model.add_transition(post, model.end, 0.1)

model.bake()

with open('exons_end_start_0.txt') as in_file:
    total = []
    for line in in_file:
        no_p_line = line.replace('P', '').replace('\n', '').lower()
        total.append(no_p_line)
converted_total = [converter_to(x, 2) for x in total]


def test(new_model):
    ok_count = 0
    wring_count = 0
    with open('exons_end_start_0.txt') as ifile:
        for line in ifile:
            end_exon = line.lower().replace('p', '').replace('\n', '')
            test_line = converter_to(end_exon, 2)

            logp, path = new_model.viterbi(test_line)
            clean_path = path[1:-1]
            print(clean_path[-51][1].name)
            if clean_path[-51][1].name == 'stop8':
                ok_count += 1
            else:
                wring_count += 1
        print(ok_count / (ok_count + wring_count))


test(model)

model.fit(converted_total,
          verbose=True,
          transition_pseudocount=1,
          emission_pseudocount=1)

test(model)

with open('partial_model_coding_to_stop_model2.json', 'w') as out:
    out.write(model.to_json())
