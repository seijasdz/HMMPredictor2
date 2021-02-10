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

matrixAcceptor0 = numpy.array(matrix_from_exa('new_acceptor1.exa'))
acceptor0_data = classify(matrixAcceptor0, 2)

model = HiddenMarkovModel('intron_acceptor')

intron = State(DiscreteDistribution(calculator.intron_calculator('cuts_intron.txt').p), name='in')
acceptor0_states = sequence_state_factory(acceptor0_data, 'acceptor0')
post = State(DiscreteDistribution(equal_distribution), name='post')

model.add_state(intron)
add_sequence(model, acceptor0_states)
model.add_state(post)

model.add_transition(model.start, intron, 1)
model.add_transition(intron, intron, 0.9)
model.add_transition(intron, acceptor0_states[0], 0.1)
model.add_transition(acceptor0_states[-1], post, 1)
model.add_transition(post, post, 0.5)
model.add_transition(post, model.end, 0.5)

model.bake()
test_l = 'GTAACACTGAATACTCAGGAACAATTAATGGATGGTAACATATGAGGAATATCTAGGAGGCACACCCTCTCTGGCATCTATGATGGGCCAAAAACCCGCATTCGCTTGGCCACAGTATGTGAAATATAACCCAGCTTAGACACAGGGTGCGGCAGCTGTCATGTTTCTCTGTGTGTGCCGAGTGTCATGTCTGCACCGTACAGGGATAGCTGAGTCTTCATCCTCCTCAGCTCCTATCTGTCCAGTGCAATGAACAGCAGCTGCTCTCTTCCTCTCTGGTTCCCATGGCAGCCATGCTCTGTTGCAGAGAGAACAGGATTGCATGTTCCCTCTTAATGGGAACGTCCATTTTGCTTTCTGGGACCACTCTCTTAATGCCGCCTGTCAAAACCAGCTAGGACTCCCTGGGGTCCAATCCCTCTGTGTTTAATCTTCTGTCATCTCTGTCCCACCTGGCTCATCAGGGAGATGCAGAAGGCTGAAGAAAAGGAAGTCCCTGAGGACTCACTGGAGGAATGTGCCATCACTTGTTCAAATAGCCATGGCCCTTATGACTCCAACCATGACTCCAACC'
converted = converter_to(test_l.lower().replace(' ', '').replace('p', ''))

#logp, path = model.viterbi(converted)
#print(logp, [x[1].name + str(i) for i, x in enumerate(path)])

with open('new_intron_acceptor.txt') as in_file:
    total = []
    for line in in_file:
        no_p_line = line.replace('P', '').replace('\n', '').lower()[-110:]
        print(len(no_p_line))
        if len(no_p_line) >= 110:
            total.append(no_p_line)
converted_total = [converter_to(x, 2) for x in total]

#print(converted_total)

def test(model):
    with open('new_intron_acceptor.txt') as test_file:
        cont_ok = 0
        cont_not_ok = 0
        for line in test_file:
                test_line = line.lower().replace('\n', '')
                logp, path = model.viterbi(converter_to(test_line, 2))
                if (path[-40][1].name == 'acceptor015'):
                    cont_ok += 1
                else:
                    cont_not_ok += 1
        print(cont_ok / (cont_not_ok + cont_ok))

test(model)

model.fit(converted_total,
          transition_pseudocount=1,
          emission_pseudocount=1,
          verbose=True)


test(model)


with open('partial_model_intron_acceptor_model.json', 'w') as out:
    out.write(model.to_json())